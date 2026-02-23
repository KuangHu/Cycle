"""Cluster transposon systems by shared protein content (protein families)."""

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse

from .config import (
    DEFAULT_COV_MODE,
    DEFAULT_COVERAGE,
    DEFAULT_FLANKING_EDIT_THRESHOLD,
    DEFAULT_JACCARD_SIM_THRESHOLD,
    DEFAULT_LOUVAIN_RESOLUTION,
    DEFAULT_MIN_SEQ_ID,
    DEFAULT_MMSEQS_THREADS,
)
from .variant_analyzer import VariantAnalyzer

log = logging.getLogger(__name__)


class SystemClusterer:
    """Cluster IS element systems by shared protein families.

    Pipeline:
      1. Extract ORF protein sequences from guide JSONs
      2. MMseqs2 clustering into protein families
      3. Build transposon x family incidence matrix
      4. Jaccard distance between transposons
      5. Louvain community detection on similarity graph
      6. Within-cluster variant analysis (flanking / guide)
    """

    def __init__(
        self,
        output_dir: str | Path,
        min_seq_id: float = DEFAULT_MIN_SEQ_ID,
        coverage: float = DEFAULT_COVERAGE,
        cov_mode: int = DEFAULT_COV_MODE,
        louvain_resolution: float = DEFAULT_LOUVAIN_RESOLUTION,
        jaccard_sim_threshold: float = DEFAULT_JACCARD_SIM_THRESHOLD,
        flanking_edit_threshold: int = DEFAULT_FLANKING_EDIT_THRESHOLD,
        mmseqs_threads: int = DEFAULT_MMSEQS_THREADS,
    ):
        self.output_dir = Path(output_dir)
        self.min_seq_id = min_seq_id
        self.coverage = coverage
        self.cov_mode = cov_mode
        self.louvain_resolution = louvain_resolution
        self.jaccard_sim_threshold = jaccard_sim_threshold
        self.flanking_edit_threshold = flanking_edit_threshold
        self.mmseqs_threads = mmseqs_threads

        self.variant_analyzer = VariantAnalyzer(
            flanking_edit_threshold=flanking_edit_threshold
        )

        # Validate mmseqs is available
        import shutil

        if shutil.which("mmseqs") is None:
            raise RuntimeError("mmseqs not found in PATH")
        log.info("mmseqs found: %s", shutil.which("mmseqs"))

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self, formatter_dirs: list[str | Path], skip_existing: bool = True) -> Path:
        """Run the full clustering pipeline.

        Parameters
        ----------
        formatter_dirs : list of paths
            Directories containing ``*/*_is_records_guide.json`` files.
        skip_existing : bool
            If True and ``system_clusters.json`` already exists, skip.

        Returns
        -------
        Path to ``system_clusters.json``.
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        results_path = self.output_dir / "system_clusters.json"

        if skip_existing and results_path.exists():
            log.info("Results already exist at %s — skipping", results_path)
            return results_path

        # 1. Extract proteins
        fasta_path, metadata = self.extract_proteins(formatter_dirs)
        if not metadata:
            log.warning("No transposon records found — nothing to cluster")
            results_path.write_text(json.dumps([], indent=2))
            return results_path

        n_with_orfs = sum(1 for m in metadata.values() if m["orfs"])
        log.info(
            "Extracted %d transposons (%d with ORFs, %d without)",
            len(metadata), n_with_orfs, len(metadata) - n_with_orfs,
        )

        # 2. Cluster proteins
        orf_to_family = self.cluster_proteins(fasta_path)

        # 3. Build incidence matrix
        incidence, trans_ids, family_ids = self.build_incidence_matrix(
            metadata, orf_to_family
        )

        # 4+5. Build similarity graph and detect communities (sparse)
        communities = self.build_graph_and_detect(incidence, trans_ids)

        # 6. Variant analysis
        clusters = self.analyze_variants(communities, metadata)

        # Write results
        results_path.write_text(json.dumps(clusters, indent=2))
        log.info("Wrote %d clusters to %s", len(clusters), results_path)

        # Write summary TSV
        self._write_summary_tsv(clusters)

        return results_path

    # ------------------------------------------------------------------
    # Step 1: Extract proteins
    # ------------------------------------------------------------------

    def extract_proteins(
        self, formatter_dirs: list[str | Path]
    ) -> tuple[Path, dict]:
        """Collect all ORF protein sequences into a master FASTA.

        Returns
        -------
        fasta_path : Path
            Path to ``all_proteins.faa``.
        metadata : dict
            Keyed by ``is_id``, each value is a dict with sample_id,
            organism, orfs, flanking_upstream, flanking_downstream,
            guide_hits, etc.
        """
        fasta_path = self.output_dir / "all_proteins.faa"
        metadata: dict[str, dict] = {}
        fasta_lines: list[str] = []

        for fmt_dir in formatter_dirs:
            fmt_dir = Path(fmt_dir)
            json_paths = sorted(fmt_dir.glob("*/*_is_records_guide.json"))
            log.info("Found %d guide JSONs in %s", len(json_paths), fmt_dir)

            for jp in json_paths:
                try:
                    records = json.loads(jp.read_text())
                except (json.JSONDecodeError, OSError) as exc:
                    log.warning("Skipping %s: %s", jp, exc)
                    continue

                for rec in records:
                    is_id = rec["is_id"]
                    if is_id in metadata:
                        continue  # deduplicate

                    orfs = (rec.get("orf_annotation") or {}).get("orfs", [])
                    orf_entries = []
                    for i, orf in enumerate(orfs):
                        seq = orf.get("protein_sequence", "")
                        if not seq:
                            continue
                        # Strip trailing stop codon marker
                        seq = seq.rstrip("*")
                        header = f"{is_id}__{orf['start']}_{orf['end']}_{orf['strand']}"
                        fasta_lines.append(f">{header}")
                        fasta_lines.append(seq)
                        orf_entries.append({
                            "header": header,
                            "start": orf["start"],
                            "end": orf["end"],
                            "strand": orf["strand"],
                            "length_nt": orf["length_nt"],
                        })

                    flank_up = rec.get("flanking_upstream", {}).get("sequence", "")
                    flank_down = rec.get("flanking_downstream", {}).get("sequence", "")

                    # Best guide hit
                    guide_hits = rec.get("guide_hits", [])
                    best_guide = None
                    if guide_hits:
                        best_guide = max(guide_hits, key=lambda h: h["length"])

                    metadata[is_id] = {
                        "sample_id": rec.get("sample_id", ""),
                        "organism": rec.get("organism", ""),
                        "orfs": orf_entries,
                        "flanking_upstream": flank_up,
                        "flanking_downstream": flank_down,
                        "guide_hits": guide_hits,
                        "best_guide_seq": best_guide["seq_noncoding"] if best_guide else "",
                        "is_length": rec.get("is_element", {}).get("length", 0),
                    }

        fasta_path.write_text("\n".join(fasta_lines) + "\n" if fasta_lines else "")
        log.info("Wrote %d protein sequences to %s", len(fasta_lines) // 2, fasta_path)
        return fasta_path, metadata

    # ------------------------------------------------------------------
    # Step 2: Cluster proteins with MMseqs2
    # ------------------------------------------------------------------

    def cluster_proteins(self, fasta_path: Path) -> dict[str, str]:
        """Run MMseqs2 easy-cluster and return orf_header -> family mapping.

        The family ID is the representative sequence header.
        """
        cluster_dir = self.output_dir / "mmseqs_clusters"
        cluster_dir.mkdir(parents=True, exist_ok=True)
        prefix = cluster_dir / "clust"
        tsv_path = Path(f"{prefix}_cluster.tsv")

        if tsv_path.exists():
            log.info("MMseqs2 cluster TSV already exists — reusing")
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                cmd = [
                    "mmseqs", "easy-cluster",
                    str(fasta_path),
                    str(prefix),
                    tmpdir,
                    "--min-seq-id", str(self.min_seq_id),
                    "-c", str(self.coverage),
                    "--cov-mode", str(self.cov_mode),
                    "--threads", str(self.mmseqs_threads),
                ]
                log.info("Running: %s", " ".join(cmd))
                ret = subprocess.run(cmd, capture_output=True, text=True)
                if ret.returncode != 0:
                    log.error("MMseqs2 stderr:\n%s", ret.stderr)
                    raise RuntimeError(f"MMseqs2 failed (exit {ret.returncode})")

        # Parse cluster TSV: representative\tmember
        orf_to_family: dict[str, str] = {}
        for line in tsv_path.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                rep, member = parts[0], parts[1]
                orf_to_family[member] = rep
        log.info(
            "MMseqs2: %d ORFs in %d families",
            len(orf_to_family), len(set(orf_to_family.values())),
        )
        return orf_to_family

    # ------------------------------------------------------------------
    # Step 3: Build incidence matrix
    # ------------------------------------------------------------------

    def build_incidence_matrix(
        self,
        metadata: dict[str, dict],
        orf_to_family: dict[str, str],
    ) -> tuple[sparse.csr_matrix, list[str], list[str]]:
        """Build a transposon x protein_family binary incidence matrix.

        Returns (csr_matrix, transposon_ids, family_ids).
        """
        # Collect family set per transposon
        all_families: set[str] = set()
        trans_families: dict[str, set[str]] = {}

        for is_id, meta in metadata.items():
            fams = set()
            for orf in meta["orfs"]:
                fam = orf_to_family.get(orf["header"])
                if fam:
                    fams.add(fam)
                    all_families.add(fam)
            trans_families[is_id] = fams

        trans_ids = sorted(trans_families.keys())
        family_ids = sorted(all_families)
        fam_to_idx = {f: i for i, f in enumerate(family_ids)}
        trans_to_idx = {t: i for i, t in enumerate(trans_ids)}

        rows, cols = [], []
        for is_id, fams in trans_families.items():
            ti = trans_to_idx[is_id]
            for f in fams:
                rows.append(ti)
                cols.append(fam_to_idx[f])

        n_trans = len(trans_ids)
        n_fam = len(family_ids)
        data = np.ones(len(rows), dtype=np.float32)
        incidence = sparse.csr_matrix(
            (data, (rows, cols)), shape=(n_trans, n_fam)
        )

        # Save artifacts
        sparse.save_npz(self.output_dir / "incidence_matrix.npz", incidence)
        (self.output_dir / "transposon_ids.json").write_text(
            json.dumps(trans_ids, indent=2)
        )
        (self.output_dir / "family_ids.json").write_text(
            json.dumps(family_ids, indent=2)
        )
        log.info("Incidence matrix: %d transposons x %d families", n_trans, n_fam)
        return incidence, trans_ids, family_ids

    # ------------------------------------------------------------------
    # Steps 4+5: Sparse Jaccard graph + community detection
    # ------------------------------------------------------------------

    def build_graph_and_detect(
        self,
        incidence: sparse.csr_matrix,
        transposon_ids: list[str],
    ) -> list[set[str]]:
        """Build a Jaccard-similarity graph (sparse) and run Louvain.

        Uses M @ M.T to find transposon pairs sharing protein families
        without materializing a dense N×N matrix.
        """
        n = len(transposon_ids)
        if n <= 1:
            return [set(transposon_ids)]

        # Number of families per transposon (for Jaccard denominator)
        row_sizes = np.asarray(incidence.sum(axis=1)).ravel()

        # Shared family counts: sparse (n x n)
        log.info("Computing sparse M @ M.T (%d transposons)...", n)
        shared = (incidence @ incidence.T).tocoo()

        # Build graph with Jaccard similarity edges
        G = nx.Graph()
        G.add_nodes_from(range(n))

        n_edges = 0
        for i, j, v in zip(shared.row, shared.col, shared.data):
            if i >= j:
                continue  # upper triangle only
            union = row_sizes[i] + row_sizes[j] - v
            if union == 0:
                continue
            sim = v / union  # Jaccard similarity
            if sim >= self.jaccard_sim_threshold:
                G.add_edge(i, j, weight=float(sim))
                n_edges += 1

        log.info("Similarity graph: %d nodes, %d edges", n, n_edges)

        communities_idx = nx.community.louvain_communities(
            G, weight="weight", resolution=self.louvain_resolution, seed=42
        )

        # Convert index sets to is_id sets
        communities = []
        for comm in communities_idx:
            communities.append({transposon_ids[i] for i in comm})

        log.info(
            "Louvain: %d communities (sizes: %s)",
            len(communities),
            ", ".join(
                str(len(c))
                for c in sorted(communities, key=len, reverse=True)[:10]
            ),
        )
        return communities

    # ------------------------------------------------------------------
    # Step 6: Variant analysis
    # ------------------------------------------------------------------

    def analyze_variants(
        self, communities: list[set[str]], metadata: dict[str, dict]
    ) -> list[dict]:
        """Run variant analysis within each community.

        Returns list of cluster dicts with L1/L2 variant info.
        """
        clusters = []
        for i, comm in enumerate(sorted(communities, key=len, reverse=True)):
            member_ids = sorted(comm)

            # Check if this is a no-ORF singleton
            all_no_orfs = all(
                not metadata[mid]["orfs"] for mid in member_ids if mid in metadata
            )

            variant_info = self.variant_analyzer.classify_cluster(
                member_ids, metadata
            )

            cluster = {
                "cluster_id": i,
                "size": len(member_ids),
                "members": member_ids,
                "no_orfs": all_no_orfs and len(member_ids) == 1,
                "variants": variant_info,
            }
            clusters.append(cluster)

        return clusters

    # ------------------------------------------------------------------
    # Summary TSV
    # ------------------------------------------------------------------

    def _write_summary_tsv(self, clusters: list[dict]) -> None:
        """Write a flat TSV with one row per transposon."""
        tsv_path = self.output_dir / "system_clusters_summary.tsv"
        header = [
            "is_id", "sample_id", "organism", "cluster_id", "cluster_size",
            "level1_variant", "level2_variant", "is_length",
            "n_orfs", "best_guide_length", "no_orfs",
        ]
        lines = ["\t".join(header)]

        for cl in clusters:
            cid = cl["cluster_id"]
            csize = cl["size"]
            variants = cl.get("variants", {})
            l1_map = variants.get("l1_assignments", {})
            l2_map = variants.get("l2_assignments", {})

            for mid in cl["members"]:
                meta = variants.get("member_meta", {}).get(mid, {})
                sample_id = meta.get("sample_id", "")
                organism = meta.get("organism", "")
                is_length = meta.get("is_length", 0)
                n_orfs = meta.get("n_orfs", 0)
                best_gl = meta.get("best_guide_length", 0)

                l1 = l1_map.get(mid, "NA")
                l2 = l2_map.get(mid, "NA")

                lines.append("\t".join(str(x) for x in [
                    mid, sample_id, organism, cid, csize,
                    l1, l2, is_length, n_orfs, best_gl, cl["no_orfs"],
                ]))

        tsv_path.write_text("\n".join(lines) + "\n")
        log.info("Wrote summary TSV: %s (%d data rows)", tsv_path, len(lines) - 1)
