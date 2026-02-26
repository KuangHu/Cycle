"""BLAST IS elements against ISfinder and score cluster novelty."""

import csv
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from .config import (
    ANNOTATED_RECORDS_JSON,
    CLUSTER_NOVELTY_JSON,
    CLUSTER_NOVELTY_TSV,
    DEFAULT_BLAST_EVALUE,
    DEFAULT_BLAST_MAX_TARGET_SEQS,
    DEFAULT_BLAST_THREADS,
    DIVERGENT_THRESHOLD,
    MIN_QUERY_COVERAGE,
    NOVEL_THRESHOLD,
    PIDENT_CEILING,
    PIDENT_FLOOR,
    WEIGHT_DIVERGENCE,
    WEIGHT_DIVERSITY,
    WEIGHT_MOSAIC,
)

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# ISfinder header parsing (reused from cycle.preprocess.is_reference)
# ------------------------------------------------------------------

def _parse_isfinder_header(header: str) -> Optional[dict]:
    """Parse ``ISname_ISgroup_ISfamily`` into components.

    Returns dict with *name*, *group*, *family* or ``None``.
    """
    parts = header.split("_")
    if len(parts) >= 3:
        family = parts[-1]
        group = parts[-2]
        name = "_".join(parts[:-2])
        return {"name": name, "group": group, "family": family}

    m = re.match(r"^(IS\S+?)_(IS\S+?)_(IS\S+)$", header)
    if m:
        return {"name": m.group(1), "group": m.group(2), "family": m.group(3)}

    logger.warning("Cannot parse ISfinder header: %s", header)
    return None


class NoveltyAnnotator:
    """Compare IS elements to ISfinder via BLAST and score cluster novelty.

    Workflow
    --------
    1. Build BLAST nucleotide DB from ISfinder FASTA
    2. Extract IS element sequences from ``*_is_records_guide.json``
    3. Run ``blastn`` against the DB
    4. Annotate each IS record with best-hit information
    5. Score each system cluster for novelty (divergence / diversity / mosaic)
    """

    def __init__(
        self,
        output_dir: str | Path = "data/novelty_output",
        evalue: float = DEFAULT_BLAST_EVALUE,
        max_target_seqs: int = DEFAULT_BLAST_MAX_TARGET_SEQS,
        threads: int = DEFAULT_BLAST_THREADS,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.evalue = evalue
        self.max_target_seqs = max_target_seqs
        self.threads = threads

        for tool in ("blastn", "makeblastdb"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

    # ==================================================================
    # Public API
    # ==================================================================

    def run(
        self,
        formatter_dirs: list[str | Path],
        clusters_path: str | Path,
        isfinder_fasta: Optional[str | Path] = None,
        skip_existing: bool = True,
    ) -> Path:
        """Full pipeline: BLAST → annotate → score.

        Parameters
        ----------
        formatter_dirs : list
            Directories containing ``*/*_is_records_guide.json``.
        clusters_path : str | Path
            Path to ``system_clusters.json``.
        isfinder_fasta : str | Path, optional
            Path to ISfinder FASTA. Auto-detected from
            ``is_reference/ISfinder_raw.fna`` next to *formatter_dirs*
            if not provided.
        skip_existing : bool
            Skip if output files already exist.

        Returns
        -------
        Path to the output directory.
        """
        summary_tsv = self.output_dir / CLUSTER_NOVELTY_TSV
        if skip_existing and summary_tsv.exists():
            logger.info("Output already exists, skipping: %s", summary_tsv)
            return self.output_dir

        # Step 1: locate / build BLAST DB
        isfinder_fasta = self._resolve_isfinder_fasta(
            isfinder_fasta, formatter_dirs
        )
        blastdb = self._build_blastdb(isfinder_fasta)

        # Step 2: extract IS sequences
        query_fasta, record_map = self._extract_sequences(formatter_dirs)
        if not record_map:
            logger.warning("No IS records found — nothing to annotate")
            return self.output_dir

        # Step 3: BLAST
        hits = self._run_blast(query_fasta, blastdb)

        # Step 4: annotate records
        annotated = self._annotate_records(record_map, hits)

        # Step 5: load clusters and score novelty
        clusters = json.loads(Path(clusters_path).read_text())
        novelty = self._score_clusters(clusters, annotated)

        # Write outputs
        self._write_outputs(annotated, novelty)

        logger.info("Novelty annotation complete → %s", self.output_dir)
        return self.output_dir

    # ==================================================================
    # Step 1 — BLAST DB
    # ==================================================================

    def _resolve_isfinder_fasta(
        self,
        explicit: Optional[str | Path],
        formatter_dirs: list[str | Path],
    ) -> Path:
        """Find the ISfinder FASTA, checking common locations."""
        if explicit:
            p = Path(explicit)
            if p.exists():
                return p
            raise FileNotFoundError(f"ISfinder FASTA not found: {p}")

        # Auto-detect: look for is_reference/ISfinder_raw.fna relative
        # to each formatter dir (typically a sibling directory)
        for fmt_dir in formatter_dirs:
            candidate = Path(fmt_dir).parent / "is_reference" / "ISfinder_raw.fna"
            if candidate.exists():
                logger.info("Auto-detected ISfinder FASTA: %s", candidate)
                return candidate

        raise FileNotFoundError(
            "Could not locate ISfinder_raw.fna. "
            "Provide --isfinder-fasta explicitly."
        )

    def _build_blastdb(self, fasta: Path) -> Path:
        """Build a BLAST nucleotide DB (skip if cached)."""
        db_prefix = self.output_dir / "isfinder_blastdb"
        nhr = Path(f"{db_prefix}.nhr")
        if nhr.exists():
            logger.info("Using cached BLAST DB: %s", db_prefix)
            return db_prefix

        cmd = [
            "makeblastdb",
            "-in", str(fasta),
            "-dbtype", "nucl",
            "-out", str(db_prefix),
        ]
        logger.info("Building BLAST DB: %s", " ".join(cmd))
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        logger.info("BLAST DB ready: %s", db_prefix)
        return db_prefix

    # ==================================================================
    # Step 2 — extract IS sequences
    # ==================================================================

    def _extract_sequences(
        self, formatter_dirs: list[str | Path]
    ) -> tuple[Path, dict]:
        """Write all IS element sequences to a FASTA and return record map.

        Returns
        -------
        (query_fasta_path, record_map)
            *record_map* maps ``is_id`` → full record dict.
        """
        fasta_path = self.output_dir / "query_is_elements.fna"
        record_map: dict[str, dict] = {}
        fasta_lines: list[str] = []

        for fmt_dir in formatter_dirs:
            fmt_dir = Path(fmt_dir)
            json_paths = sorted(fmt_dir.glob("*/*_is_records_guide.json"))
            logger.info(
                "Found %d guide JSONs in %s", len(json_paths), fmt_dir
            )

            for jp in json_paths:
                try:
                    records = json.loads(jp.read_text())
                except (json.JSONDecodeError, OSError) as exc:
                    logger.warning("Skipping %s: %s", jp, exc)
                    continue

                for rec in records:
                    is_id = rec.get("is_id", "")
                    if not is_id or is_id in record_map:
                        continue

                    seq = (rec.get("is_element") or {}).get("sequence", "")
                    if not seq:
                        logger.debug("No sequence for %s, skipping", is_id)
                        continue

                    record_map[is_id] = rec
                    fasta_lines.append(f">{is_id}")
                    fasta_lines.append(seq)

        fasta_path.write_text(
            "\n".join(fasta_lines) + "\n" if fasta_lines else ""
        )
        logger.info(
            "Extracted %d IS sequences → %s", len(record_map), fasta_path
        )
        return fasta_path, record_map

    # ==================================================================
    # Step 3 — BLAST
    # ==================================================================

    _BLAST_OUTFMT = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen"
    )

    def _run_blast(self, query: Path, db: Path) -> dict[str, list[dict]]:
        """Run blastn and return hits grouped by query IS id.

        Returns dict mapping ``is_id`` → list of hit dicts, sorted by
        descending bitscore.
        """
        out_path = self.output_dir / "blast_results.tsv"
        cmd = [
            "blastn",
            "-query", str(query),
            "-db", str(db),
            "-outfmt", self._BLAST_OUTFMT,
            "-evalue", str(self.evalue),
            "-num_threads", str(self.threads),
            "-max_target_seqs", str(self.max_target_seqs),
            "-out", str(out_path),
        ]
        logger.info("Running BLAST: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

        return self._parse_blast(out_path)

    def _parse_blast(self, blast_tsv: Path) -> dict[str, list[dict]]:
        """Parse BLAST tabular output into grouped hits."""
        fields = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen",
        ]
        float_fields = {"pident", "evalue", "bitscore"}
        int_fields = {
            "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "qlen", "slen",
        }

        hits: dict[str, list[dict]] = {}

        if not blast_tsv.exists() or blast_tsv.stat().st_size == 0:
            logger.info("No BLAST hits")
            return hits

        with open(blast_tsv) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < len(fields):
                    continue

                row = {}
                for i, fname in enumerate(fields):
                    val = parts[i]
                    if fname in float_fields:
                        row[fname] = float(val)
                    elif fname in int_fields:
                        row[fname] = int(val)
                    else:
                        row[fname] = val

                # Compute query coverage
                aln_span = abs(row["qend"] - row["qstart"]) + 1
                row["query_coverage"] = aln_span / row["qlen"] if row["qlen"] else 0.0

                # Parse ISfinder subject header
                subj_header = row["sseqid"].split()[0]
                parsed = _parse_isfinder_header(subj_header)
                row["hit_name"] = parsed["name"] if parsed else subj_header
                row["hit_group"] = parsed["group"] if parsed else ""
                row["hit_family"] = parsed["family"] if parsed else ""

                hits.setdefault(row["qseqid"], []).append(row)

        # Sort each query's hits by bitscore descending
        for qid in hits:
            hits[qid].sort(key=lambda h: h["bitscore"], reverse=True)

        total = sum(len(v) for v in hits.values())
        logger.info(
            "Parsed %d BLAST hits for %d queries", total, len(hits)
        )
        return hits

    # ==================================================================
    # Step 4 — annotate records
    # ==================================================================

    def _annotate_records(
        self,
        record_map: dict[str, dict],
        hits: dict[str, list[dict]],
    ) -> dict[str, dict]:
        """Add ``isfinder_annotation`` to each IS record.

        Returns the updated record_map (mutated in-place).
        """
        n_with_hit = 0

        for is_id, rec in record_map.items():
            query_hits = hits.get(is_id, [])
            # Filter by minimum query coverage
            query_hits = [
                h for h in query_hits
                if h["query_coverage"] >= MIN_QUERY_COVERAGE
            ]

            if not query_hits:
                rec["isfinder_annotation"] = {
                    "best_hit_name": None,
                    "best_hit_family": None,
                    "pident": None,
                    "query_coverage": None,
                    "evalue": None,
                    "bitscore": None,
                    "all_hit_families": [],
                    "all_hit_names": [],
                    "is_multi_family": False,
                    "n_hits": 0,
                }
                continue

            best = query_hits[0]
            all_families = list(
                dict.fromkeys(h["hit_family"] for h in query_hits if h["hit_family"])
            )
            all_names = list(
                dict.fromkeys(h["hit_name"] for h in query_hits if h["hit_name"])
            )

            rec["isfinder_annotation"] = {
                "best_hit_name": best["hit_name"],
                "best_hit_family": best["hit_family"],
                "pident": best["pident"],
                "query_coverage": round(best["query_coverage"], 4),
                "evalue": best["evalue"],
                "bitscore": best["bitscore"],
                "all_hit_families": all_families,
                "all_hit_names": all_names,
                "is_multi_family": len(all_families) > 1,
                "n_hits": len(query_hits),
            }
            n_with_hit += 1

        logger.info(
            "Annotated %d/%d IS records with ISfinder hits",
            n_with_hit, len(record_map),
        )
        return record_map

    # ==================================================================
    # Step 5 — cluster novelty scoring
    # ==================================================================

    def _score_clusters(
        self,
        clusters: list[dict],
        record_map: dict[str, dict],
    ) -> list[dict]:
        """Compute novelty scores for each system cluster."""
        results = []

        for cl in clusters:
            cluster_id = cl["cluster_id"]
            members = cl.get("members", [])
            variants = cl.get("variants") or {}

            # Gather annotations for members present in record_map
            member_annotations = {}
            for mid in members:
                rec = record_map.get(mid)
                if rec and "isfinder_annotation" in rec:
                    member_annotations[mid] = rec["isfinder_annotation"]

            # Basic stats
            n_with_match = sum(
                1 for a in member_annotations.values()
                if a.get("best_hit_name") is not None
            )
            n_without_match = len(member_annotations) - n_with_match

            pidents = [
                a["pident"]
                for a in member_annotations.values()
                if a.get("pident") is not None
            ]
            mean_pident = sum(pidents) / len(pidents) if pidents else 0.0
            min_pident = min(pidents) if pidents else 0.0

            # Collect all ISfinder families and names across members
            all_families: set[str] = set()
            all_names: set[str] = set()
            for a in member_annotations.values():
                all_families.update(a.get("all_hit_families", []))
                all_names.update(a.get("all_hit_names", []))

            # Dominant family/name (most common best hit)
            fam_counts: dict[str, int] = {}
            name_counts: dict[str, int] = {}
            for a in member_annotations.values():
                f = a.get("best_hit_family")
                n = a.get("best_hit_name")
                if f:
                    fam_counts[f] = fam_counts.get(f, 0) + 1
                if n:
                    name_counts[n] = name_counts.get(n, 0) + 1

            dominant_family = (
                max(fam_counts, key=fam_counts.get) if fam_counts else ""
            )
            dominant_name = (
                max(name_counts, key=name_counts.get) if name_counts else ""
            )

            # --- Divergence score (weight 0.5) ---
            divergence = self._divergence_score(
                n_with_match, n_without_match, mean_pident
            )

            # --- Diversity score (weight 0.3) ---
            n_l1 = variants.get("n_l1_groups", 1)
            n_l2 = variants.get("n_l2_groups", 1)
            diversity = self._diversity_score(n_l1, n_l2, len(members))

            # --- Mosaic score (weight 0.2) ---
            mosaic = self._mosaic_score(all_families, all_names)

            # --- Composite ---
            novelty_score = (
                WEIGHT_DIVERGENCE * divergence
                + WEIGHT_DIVERSITY * diversity
                + WEIGHT_MOSAIC * mosaic
            )
            novelty_score = round(novelty_score, 4)

            if novelty_score >= NOVEL_THRESHOLD:
                novelty_class = "novel"
            elif novelty_score >= DIVERGENT_THRESHOLD:
                novelty_class = "divergent"
            else:
                novelty_class = "known"

            # Check circle evidence across members
            has_circle = any(
                (record_map.get(mid, {}).get("circle_evidence") or {})
                .get("n_tail_head_reads", 0) > 0
                for mid in members
            )

            row = {
                "cluster_id": cluster_id,
                "cluster_size": len(members),
                "n_with_match": n_with_match,
                "n_without_match": n_without_match,
                "mean_pident": round(mean_pident, 2),
                "min_pident": round(min_pident, 2),
                "dominant_isfinder_family": dominant_family,
                "dominant_isfinder_name": dominant_name,
                "n_distinct_families": len(all_families),
                "n_distinct_names": len(all_names),
                "n_l1_groups": n_l1,
                "n_l2_groups": n_l2,
                "divergence_score": round(divergence, 4),
                "diversity_score": round(diversity, 4),
                "mosaic_score": round(mosaic, 4),
                "novelty_score": novelty_score,
                "novelty_class": novelty_class,
                "has_circle_evidence": has_circle,
                "member_annotations": {
                    mid: member_annotations.get(mid)
                    for mid in members
                    if mid in member_annotations
                },
            }
            results.append(row)

        # Sort by novelty score descending
        results.sort(key=lambda r: r["novelty_score"], reverse=True)

        n_novel = sum(1 for r in results if r["novelty_class"] == "novel")
        n_divergent = sum(
            1 for r in results if r["novelty_class"] == "divergent"
        )
        logger.info(
            "Scored %d clusters: %d novel, %d divergent, %d known",
            len(results), n_novel, n_divergent,
            len(results) - n_novel - n_divergent,
        )
        return results

    # ------------------------------------------------------------------
    # Sub-scores
    # ------------------------------------------------------------------

    @staticmethod
    def _divergence_score(
        n_with: int, n_without: int, mean_pident: float
    ) -> float:
        """Score 0–1: how different from ISfinder.

        No matches → 1.0.  High %identity → 0.0.
        """
        total = n_with + n_without
        if total == 0:
            return 1.0

        # Fraction without any BLAST hit
        frac_no_hit = n_without / total

        if n_with == 0:
            return 1.0

        # Linear scale between PIDENT_FLOOR and PIDENT_CEILING
        if mean_pident >= PIDENT_CEILING:
            pident_component = 0.0
        elif mean_pident <= PIDENT_FLOOR:
            pident_component = 1.0
        else:
            pident_component = (
                (PIDENT_CEILING - mean_pident)
                / (PIDENT_CEILING - PIDENT_FLOOR)
            )

        # Blend: weight the no-hit fraction and pident-based divergence
        score = 0.5 * frac_no_hit + 0.5 * pident_component
        return min(1.0, max(0.0, score))

    @staticmethod
    def _diversity_score(
        n_l1: int, n_l2: int, cluster_size: int
    ) -> float:
        """Score 0–1: active diversification within the cluster.

        Based on L1 (flanking) and L2 (guide) variant groups relative
        to cluster size.
        """
        if cluster_size <= 1:
            return 0.0

        # Ratio of variant groups to cluster size (capped at 1.0)
        l1_ratio = min(1.0, (n_l1 - 1) / (cluster_size - 1)) if cluster_size > 1 else 0.0
        l2_ratio = min(1.0, (n_l2 - 1) / (cluster_size - 1)) if cluster_size > 1 else 0.0

        # L2 diversity (guide-level) is more informative than L1
        score = 0.4 * l1_ratio + 0.6 * l2_ratio
        return min(1.0, max(0.0, score))

    @staticmethod
    def _mosaic_score(families: set[str], names: set[str]) -> float:
        """Score 0–1: members match different ISfinder entries.

        Multiple families → 1.0.
        Multiple names in same family → 0.3–0.7 (scaled by name count).
        Single name → 0.0.
        """
        n_fam = len(families)
        n_name = len(names)

        if n_fam == 0 and n_name == 0:
            # No hits at all — ambiguous, slight positive signal
            return 0.5

        if n_fam > 1:
            return 1.0

        if n_name > 1:
            # Same family but different IS entries → moderate mosaic
            return min(0.7, 0.3 + 0.1 * (n_name - 1))

        return 0.0

    # ==================================================================
    # Output
    # ==================================================================

    def _write_outputs(
        self,
        record_map: dict[str, dict],
        novelty: list[dict],
    ) -> None:
        """Write annotated records (JSON), cluster novelty (TSV + JSON)."""

        # --- annotated_is_records.json ---
        records_path = self.output_dir / ANNOTATED_RECORDS_JSON
        records_list = list(record_map.values())
        records_path.write_text(json.dumps(records_list, indent=2))
        logger.info(
            "Wrote %d annotated records → %s", len(records_list), records_path
        )

        # --- cluster_novelty.json (includes per-member annotations) ---
        json_path = self.output_dir / CLUSTER_NOVELTY_JSON
        json_path.write_text(json.dumps(novelty, indent=2))
        logger.info("Wrote cluster novelty JSON → %s", json_path)

        # --- cluster_novelty_summary.tsv ---
        tsv_path = self.output_dir / CLUSTER_NOVELTY_TSV
        tsv_fields = [
            "cluster_id", "cluster_size",
            "n_with_match", "n_without_match",
            "mean_pident", "min_pident",
            "dominant_isfinder_family", "dominant_isfinder_name",
            "n_distinct_families", "n_distinct_names",
            "n_l1_groups", "n_l2_groups",
            "divergence_score", "diversity_score", "mosaic_score",
            "novelty_score", "novelty_class",
            "has_circle_evidence",
        ]
        with open(tsv_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=tsv_fields, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            for row in novelty:
                writer.writerow(row)

        logger.info(
            "Wrote %d rows → %s", len(novelty), tsv_path
        )
