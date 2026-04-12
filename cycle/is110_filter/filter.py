"""
IS110 Element Filter

Identifies IS elements containing IS110-family transposase (DEDD + Tnp20 domains).
Splits results into two groups:
  - with_circle_evidence: tail-head junction reads > 0 (or partial circle detected)
  - without_circle_evidence: IS110 transposase but no circle evidence

Workflow:
  1. Run hmmsearch with DEDD.hmm and Tnp20.hmm against a protein FASTA
  2. Parse hits at protein level, require both domains per transposon
  3. Scan *_is_records_guide.json files, split IS110 records by circle evidence
  4. Load partial circle data for additional circle evidence
  5. Export JSON + TSV summaries and generate PNG + GBK visualizations
"""

import csv
import json
import logging
import os
import subprocess
from collections import defaultdict
from glob import glob
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

DEFAULT_DEDD_HMM = "/home/kuangh/scripts/IS110/hmm/DEDD.hmm"
DEFAULT_TNP20_HMM = "/home/kuangh/scripts/IS110/hmm/Tnp20.hmm"


class IS110Filter:
    """Filter IS elements for IS110 transposase, split by circle evidence."""

    def __init__(
        self,
        dedd_hmm: str = DEFAULT_DEDD_HMM,
        tnp20_hmm: str = DEFAULT_TNP20_HMM,
        evalue: float = 1e-5,
        max_orf_gap: int = 300,
    ):
        self.dedd_hmm = dedd_hmm
        self.tnp20_hmm = tnp20_hmm
        self.evalue = evalue
        self.max_orf_gap = max_orf_gap  # max bp gap between ORFs for case 2

    # ------------------------------------------------------------------
    # Step 1: hmmsearch
    # ------------------------------------------------------------------

    def hmmsearch(
        self,
        protein_fasta: str,
        output_tbl: str,
        hmm_path: str,
        cpus: int = 8,
    ) -> str:
        """Run hmmsearch and return path to tblout file."""
        cmd = [
            "hmmsearch",
            "--cpu", str(cpus),
            "--tblout", output_tbl,
            "-E", str(self.evalue),
            hmm_path,
            protein_fasta,
        ]
        logger.info("Running: %s", " ".join(cmd))
        subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
        return output_tbl

    def run_hmmsearch(
        self,
        protein_fasta: str,
        output_dir: str,
        cpus: int = 8,
    ) -> Tuple[Set[str], Dict[str, Dict[str, List[Tuple[str, int, int, str]]]]]:
        """Run both DEDD and Tnp20 hmmsearch, return IS110 transposon IDs and domain map.

        IS110 requires BOTH DEDD and Tnp20 domains:
          Case 1: single ORF hit by both HMMs
          Case 2: two adjacent ORFs (gap <= max_orf_gap), one with DEDD, one with Tnp20

        Returns:
            (is110_ids, trans_domains) where trans_domains maps
            transposon_id -> {"DEDD": [...], "Tnp20": [...]} with ORF coords.
        """
        os.makedirs(output_dir, exist_ok=True)

        # Run hmmsearch for each domain, collect per-protein hits
        protein_hits: Dict[str, Dict[str, Tuple[str, int, int, str]]] = {}
        for name, hmm in [("DEDD", self.dedd_hmm), ("Tnp20", self.tnp20_hmm)]:
            tbl = os.path.join(output_dir, f"{name}_hits.tbl")
            self.hmmsearch(protein_fasta, tbl, hmm, cpus=cpus)
            hits = self.parse_tblout(tbl)
            protein_hits[name] = hits
            logger.info(
                "%s: %d protein hits → %d unique transposons",
                name, len(hits),
                len({v[0] for v in hits.values()}),
            )

        return self._identify_is110(protein_hits)

    def _identify_is110(
        self, protein_hits: Dict[str, Dict[str, Tuple[str, int, int, str]]],
    ) -> Tuple[Set[str], Dict[str, Dict[str, List[Tuple[str, int, int, str]]]]]:
        """Apply two-case logic to identify IS110 transposons from HMM hits."""
        # Group proteins by transposon and domain
        trans_domains: Dict[str, Dict[str, List[Tuple[str, int, int, str]]]] = defaultdict(
            lambda: {"DEDD": [], "Tnp20": []}
        )
        for domain, hits in protein_hits.items():
            for prot_id, (trans_id, start, end, strand) in hits.items():
                trans_domains[trans_id][domain].append((prot_id, start, end, strand))

        is110_ids: Set[str] = set()
        n_case1 = 0
        n_case2 = 0

        for trans_id, domains in trans_domains.items():
            dedd_prots = {p[0] for p in domains["DEDD"]}
            tnp20_prots = {p[0] for p in domains["Tnp20"]}

            # Case 1: same ORF has both domains
            if dedd_prots & tnp20_prots:
                is110_ids.add(trans_id)
                n_case1 += 1
                continue

            # Case 2: two adjacent ORFs, one DEDD, one Tnp20
            if not dedd_prots or not tnp20_prots:
                continue

            if self._has_adjacent_pair(domains["DEDD"], domains["Tnp20"]):
                is110_ids.add(trans_id)
                n_case2 += 1

        logger.info(
            "IS110: %d transposons (case1=%d same-ORF, case2=%d adjacent-ORFs, "
            "rejected %d single-domain)",
            len(is110_ids), n_case1, n_case2,
            len(trans_domains) - len(is110_ids),
        )
        return is110_ids, dict(trans_domains)

    def _has_adjacent_pair(
        self,
        dedd_orfs: List[Tuple[str, int, int, str]],
        tnp20_orfs: List[Tuple[str, int, int, str]],
    ) -> bool:
        """Check if any DEDD ORF is adjacent to any Tnp20 ORF (gap <= max_orf_gap)."""
        for _, ds, de, _ in dedd_orfs:
            for _, ts, te, _ in tnp20_orfs:
                # Gap between the two ORFs (0 if overlapping)
                if de <= ts:
                    gap = ts - de
                elif te <= ds:
                    gap = ds - te
                else:
                    gap = 0  # overlapping
                if gap <= self.max_orf_gap:
                    return True
        return False

    # ------------------------------------------------------------------
    # Step 2: parse hmmsearch output
    # ------------------------------------------------------------------

    def parse_tblout(
        self, tbl_path: str,
    ) -> Dict[str, Tuple[str, int, int, str]]:
        """Parse hmmsearch tblout, return per-protein hit info.

        Header format: {transposon_id}__{start}_{end}_{strand}

        Returns:
            dict mapping protein_id -> (transposon_id, start, end, strand)
        """
        hits: Dict[str, Tuple[str, int, int, str]] = {}
        with open(tbl_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                protein_id = line.split()[0]
                if protein_id in hits:
                    continue  # keep first (best E-value)
                transposon_id, coords = protein_id.rsplit("__", 1)
                parts = coords.split("_")
                start, end, strand = int(parts[0]), int(parts[1]), parts[2]
                hits[protein_id] = (transposon_id, start, end, strand)
        return hits

    # ------------------------------------------------------------------
    # Step 3: load partial circle data
    # ------------------------------------------------------------------

    @staticmethod
    def load_partial_circles(partial_circle_dirs: List[str]) -> Dict[str, List[Dict]]:
        """Load partial circle calls, keyed by is_id.

        Args:
            partial_circle_dirs: list of directories containing
                *_partial_circle_summary.json files (searched recursively).

        Returns:
            dict mapping is_id -> list of partial circle call dicts.
        """
        pc_by_is: Dict[str, List[Dict]] = {}
        for pc_dir in partial_circle_dirs:
            if not os.path.isdir(pc_dir):
                continue
            for pcf in glob(os.path.join(pc_dir, "**/*_partial_circle_summary.json"),
                            recursive=True):
                with open(pcf) as f:
                    data = json.load(f)
                if isinstance(data, list):
                    for call in data:
                        is_id = call.get("is_id", "")
                        if is_id:
                            pc_by_is.setdefault(is_id, []).append(call)
        logger.info("Loaded partial circle data for %d IS elements", len(pc_by_is))
        return pc_by_is

    # ------------------------------------------------------------------
    # Step 4: filter and split records
    # ------------------------------------------------------------------

    def filter_records(
        self,
        formatter_dir: str,
        is110_ids: Set[str],
        trans_domains: Optional[Dict] = None,
        partial_circle_ids: Optional[Set[str]] = None,
    ) -> Tuple[List[Dict], List[Dict]]:
        """Scan IS records, split IS110 by circle evidence.

        Scans guide JSONs for records with circle evidence, then scans
        annotated JSONs to find IS110 records without circle evidence
        (which were excluded from guide JSONs by the TH filter).

        Circle evidence = tail-head reads > 0 OR partial circle detected.

        Returns:
            (with_circle, without_circle) — two lists of record dicts.
        """
        if partial_circle_ids is None:
            partial_circle_ids = set()

        with_circle = []
        without_circle = []
        seen_ids: Set[str] = set()

        # Pass 1: guide JSONs — these have circle evidence (TH > 0)
        guide_files = sorted(glob(os.path.join(formatter_dir, "*", "*_is_records_guide.json")))
        for jf in guide_files:
            with open(jf) as f:
                records = json.load(f)
            for rec in records:
                is_id = rec.get("is_id", "")
                if is_id not in is110_ids:
                    continue
                if trans_domains:
                    self._annotate_orf_domains(rec, trans_domains)
                ce = rec.get("circle_evidence", {})
                n_th = ce.get("n_tail_head_reads", 0)
                has_full = n_th > 0
                has_partial = is_id in partial_circle_ids
                if has_full or has_partial:
                    rec["_circle_type"] = "full" if has_full else "partial"
                    with_circle.append(rec)
                else:
                    rec["_circle_type"] = "none"
                    without_circle.append(rec)
                seen_ids.add(is_id)

        # Pass 2: annotated JSONs — pick up IS110 records excluded by TH filter
        ann_files = sorted(glob(os.path.join(formatter_dir, "*", "*_is_records_annotated.json")))
        for jf in ann_files:
            with open(jf) as f:
                records = json.load(f)
            for rec in records:
                is_id = rec.get("is_id", "")
                if is_id not in is110_ids or is_id in seen_ids:
                    continue
                if trans_domains:
                    self._annotate_orf_domains(rec, trans_domains)
                # Check if partial circle gives it circle evidence
                if is_id in partial_circle_ids:
                    rec["_circle_type"] = "partial"
                    with_circle.append(rec)
                else:
                    rec["_circle_type"] = "none"
                    without_circle.append(rec)
                seen_ids.add(is_id)

        logger.info(
            "Scanned %s: %d IS110 total, %d with circle evidence, %d without",
            os.path.basename(formatter_dir), len(with_circle) + len(without_circle),
            len(with_circle), len(without_circle),
        )
        return with_circle, without_circle

    @staticmethod
    def _annotate_orf_domains(rec: Dict, trans_domains: Dict) -> None:
        """Add 'domains' list to each ORF that has DEDD/Tnp20 hits."""
        is_id = rec.get("is_id", "")
        domains_info = trans_domains.get(is_id)
        if not domains_info:
            return

        orfs = rec.get("orf_annotation", {}).get("orfs", [])
        coord_domains: Dict[Tuple[int, int], List[str]] = defaultdict(list)
        for domain_name in ("DEDD", "Tnp20"):
            for _prot_id, start, end, _strand in domains_info.get(domain_name, []):
                coord_domains[(start, end)].append(domain_name)

        for orf in orfs:
            orf_key = (orf["start"], orf["end"])
            if orf_key in coord_domains:
                orf["domains"] = sorted(set(coord_domains[orf_key]))

    # ------------------------------------------------------------------
    # Step 5: export
    # ------------------------------------------------------------------

    def export_results(
        self,
        with_circle: List[Dict],
        without_circle: List[Dict],
        output_dir: str,
    ):
        """Write split results as JSON and summary TSV.

        Outputs:
            {output_dir}/with_circle_evidence/
                is110_records.json + is110_summary.tsv
            {output_dir}/without_circle_evidence/
                is110_records.json + is110_summary.tsv
        """
        for label, records in [
            ("with_circle_evidence", with_circle),
            ("without_circle_evidence", without_circle),
        ]:
            sub_dir = os.path.join(output_dir, label)
            os.makedirs(sub_dir, exist_ok=True)
            self._write_json_and_tsv(records, sub_dir)

    @staticmethod
    def _write_json_and_tsv(records: List[Dict], output_dir: str):
        """Write records as JSON and summary TSV to output_dir."""
        json_path = os.path.join(output_dir, "is110_records.json")
        with open(json_path, "w") as f:
            json.dump(records, f, indent=2)
        logger.info("Wrote %d records to %s", len(records), json_path)

        tsv_path = os.path.join(output_dir, "is110_summary.tsv")
        fields = [
            "is_id", "sample_id", "is_length", "n_orfs", "circle_type",
            "n_tail_head_reads", "n_genome_head_reads", "n_tail_genome_reads",
            "n_guide_hits", "best_guide_length",
        ]
        with open(tsv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for rec in records:
                ce = rec.get("circle_evidence", {})
                gs = rec.get("guide_summary", {})
                writer.writerow({
                    "is_id": rec.get("is_id", ""),
                    "sample_id": rec.get("sample_id", ""),
                    "is_length": rec.get("is_element", {}).get("length", 0),
                    "n_orfs": (rec.get("orf_annotation") or {}).get("num_orfs", 0),
                    "circle_type": rec.get("_circle_type", ""),
                    "n_tail_head_reads": ce.get("n_tail_head_reads", 0),
                    "n_genome_head_reads": ce.get("n_genome_head_reads", 0),
                    "n_tail_genome_reads": ce.get("n_tail_genome_reads", 0),
                    "n_guide_hits": gs.get("n_hits", 0),
                    "best_guide_length": gs.get("best_length", 0),
                })
        logger.info("Wrote summary to %s", tsv_path)

    # ------------------------------------------------------------------
    # Step 6: visualize
    # ------------------------------------------------------------------

    @staticmethod
    def visualize(
        records: List[Dict],
        output_dir: str,
        pc_by_is: Optional[Dict[str, List[Dict]]] = None,
        dpi: int = 150,
    ) -> Tuple[int, int]:
        """Generate PNG + GBK for each record.

        Args:
            records: list of IS110 record dicts.
            output_dir: directory for PNG and GBK files.
            pc_by_is: partial circle data keyed by is_id.
            dpi: PNG resolution.

        Returns:
            (n_done, n_errors)
        """
        from cycle.visualizer.visualizer import ISElementVisualizer
        from cycle.visualizer.genbank import ISElementGenBank

        os.makedirs(output_dir, exist_ok=True)
        vis = ISElementVisualizer()
        gb = ISElementGenBank()

        if pc_by_is is None:
            pc_by_is = {}

        n_done = 0
        n_err = 0

        for rec in records:
            is_id = rec["is_id"]
            sample_id = rec["sample_id"]
            tag = f"{sample_id}__{is_id[:8]}"
            png_path = os.path.join(output_dir, f"{tag}.png")
            gbk_path = os.path.join(output_dir, f"{tag}.gbk")
            alignments = rec.get("guide_hits", [])
            circle_info = rec.get("circle_evidence")
            partial_circles = pc_by_is.get(is_id)

            try:
                vis.save_element_png(
                    rec, alignments, png_path, dpi=dpi,
                    circle_info=circle_info, partial_circles=partial_circles,
                )
                gb.save_genbank(rec, alignments, gbk_path, circle_info=circle_info)
                n_done += 1
            except Exception as e:
                logger.error("Failed %s: %s", tag, e)
                n_err += 1

            if n_done % 100 == 0 and n_done > 0:
                logger.info("  %d/%d done", n_done, len(records))

        logger.info("Visualized: %d PNG+GBK, %d errors", n_done, n_err)
        return n_done, n_err

    # ------------------------------------------------------------------
    # All-in-one
    # ------------------------------------------------------------------

    def run(
        self,
        protein_fasta: str,
        formatter_dirs: List[str],
        output_dir: str,
        cpus: int = 8,
        partial_circle_dirs: Optional[List[str]] = None,
        visualize: bool = True,
        dpi: int = 150,
    ) -> Tuple[List[Dict], List[Dict]]:
        """Full pipeline: hmmsearch → filter → export → visualize.

        Args:
            protein_fasta: path to all_proteins.faa.
            formatter_dirs: list of formatter output directories.
            output_dir: directory for output files.
            cpus: threads for hmmsearch.
            partial_circle_dirs: directories with partial circle summaries.
            visualize: whether to generate PNG + GBK.
            dpi: PNG resolution.

        Returns:
            (with_circle, without_circle) record lists.
        """
        is110_ids, trans_domains = self.run_hmmsearch(protein_fasta, output_dir, cpus=cpus)

        # Load partial circle data
        pc_by_is: Dict[str, List[Dict]] = {}
        partial_circle_ids: Set[str] = set()
        if partial_circle_dirs:
            pc_by_is = self.load_partial_circles(partial_circle_dirs)
            partial_circle_ids = set(pc_by_is.keys())

        # Filter and split across all formatter dirs
        all_with = []
        all_without = []
        for fmt_dir in formatter_dirs:
            w, wo = self.filter_records(
                fmt_dir, is110_ids,
                trans_domains=trans_domains,
                partial_circle_ids=partial_circle_ids,
            )
            all_with.extend(w)
            all_without.extend(wo)

        logger.info(
            "Total: %d with circle evidence, %d without",
            len(all_with), len(all_without),
        )

        # Export
        self.export_results(all_with, all_without, output_dir)

        # Visualize
        if visualize:
            wc_dir = os.path.join(output_dir, "with_circle_evidence")
            wo_dir = os.path.join(output_dir, "without_circle_evidence")
            self.visualize(all_with, wc_dir, pc_by_is=pc_by_is, dpi=dpi)
            self.visualize(all_without, wo_dir, pc_by_is=pc_by_is, dpi=dpi)

        return all_with, all_without
