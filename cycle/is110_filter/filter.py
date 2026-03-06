"""
IS110 Circular Element Filter

Identifies IS elements that:
  1. Have tail-head junction reads (circle_evidence.n_tail_head_reads > 0)
  2. Contain IS110-family transposase — requires BOTH DEDD and Tnp20 domains:
     a. Single ORF with both DEDD and Tnp20 HMM hits, OR
     b. Two adjacent ORFs (within max_orf_gap bp), one with DEDD and one with Tnp20

Workflow:
  1. Run hmmsearch with DEDD.hmm and Tnp20.hmm against a protein FASTA
  2. Parse hits at protein level, require both domains per transposon
  3. Scan *_is_records_guide.json files, emit records matching both criteria
"""

import csv
import json
import logging
import os
import subprocess
import tempfile
from collections import defaultdict
from glob import glob
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

DEFAULT_DEDD_HMM = "/home/kuangh/scripts/IS110/hmm/DEDD.hmm"
DEFAULT_TNP20_HMM = "/home/kuangh/scripts/IS110/hmm/Tnp20.hmm"


class IS110Filter:
    """Filter IS elements for IS110 transposase + tail-head circular junctions."""

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
        # protein_hits[domain] = {protein_id: (transposon_id, start, end, strand)}
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
    # Step 3: filter records
    # ------------------------------------------------------------------

    def filter_records(
        self,
        formatter_dir: str,
        is110_ids: Set[str],
        min_tail_head: int = 1,
        trans_domains: Optional[Dict] = None,
    ) -> List[Dict]:
        """Scan guide JSONs, return records with IS110 protein AND tail-head reads.

        Args:
            formatter_dir: directory with sample subdirs containing *_is_records_guide.json.
            is110_ids: set of transposon IDs that contain IS110 protein.
            min_tail_head: minimum n_tail_head_reads to pass filter.
            trans_domains: if provided, annotate each ORF with its domain hits.

        Returns:
            List of matching record dicts (with domain annotations if trans_domains given).
        """
        json_files = sorted(glob(os.path.join(formatter_dir, "*", "*_is_records_guide.json")))
        if not json_files:
            logger.warning("No guide JSON files found under %s", formatter_dir)
            return []

        matched = []
        total_scanned = 0
        total_is110 = 0
        total_circular = 0

        for jf in json_files:
            with open(jf) as f:
                records = json.load(f)

            for rec in records:
                total_scanned += 1
                is_id = rec.get("is_id", "")
                has_is110 = is_id in is110_ids
                ce = rec.get("circle_evidence", {})
                n_th = ce.get("n_tail_head_reads", 0)
                has_circle = n_th >= min_tail_head

                if has_is110:
                    total_is110 += 1
                if has_circle:
                    total_circular += 1
                if has_is110 and has_circle:
                    if trans_domains:
                        self._annotate_orf_domains(rec, trans_domains)
                    matched.append(rec)

        logger.info(
            "Scanned %d records: %d IS110, %d circular (TH>=%d), %d both",
            total_scanned, total_is110, total_circular, min_tail_head, len(matched),
        )
        return matched

    @staticmethod
    def _annotate_orf_domains(rec: Dict, trans_domains: Dict) -> None:
        """Add 'domains' list to each ORF that has DEDD/Tnp20 hits.

        Matches protein hits to ORFs by comparing (start, end) coordinates.
        The protein header encodes ORF coords as {is_id}__{start}_{end}_{strand}.
        """
        is_id = rec.get("is_id", "")
        domains_info = trans_domains.get(is_id)
        if not domains_info:
            return

        orfs = rec.get("orf_annotation", {}).get("orfs", [])
        # Build lookup: (start, end) -> set of domain names
        coord_domains: Dict[Tuple[int, int], List[str]] = defaultdict(list)
        for domain_name in ("DEDD", "Tnp20"):
            for _prot_id, start, end, _strand in domains_info.get(domain_name, []):
                coord_domains[(start, end)].append(domain_name)

        for orf in orfs:
            orf_key = (orf["start"], orf["end"])
            if orf_key in coord_domains:
                orf["domains"] = sorted(set(coord_domains[orf_key]))

    # ------------------------------------------------------------------
    # Step 4: export
    # ------------------------------------------------------------------

    def export_results(
        self,
        records: List[Dict],
        output_dir: str,
    ):
        """Write filtered records as JSON and summary TSV.

        Outputs:
            {output_dir}/is110_circular_records.json  — full records
            {output_dir}/is110_circular_summary.tsv   — one-line-per-record summary
        """
        os.makedirs(output_dir, exist_ok=True)

        json_path = os.path.join(output_dir, "is110_circular_records.json")
        with open(json_path, "w") as f:
            json.dump(records, f, indent=2)
        logger.info("Wrote %d records to %s", len(records), json_path)

        tsv_path = os.path.join(output_dir, "is110_circular_summary.tsv")
        fields = [
            "is_id", "sample_id", "is_length", "n_orfs",
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
                    "n_orfs": rec.get("orf_annotation", {}).get("num_orfs", 0),
                    "n_tail_head_reads": ce.get("n_tail_head_reads", 0),
                    "n_genome_head_reads": ce.get("n_genome_head_reads", 0),
                    "n_tail_genome_reads": ce.get("n_tail_genome_reads", 0),
                    "n_guide_hits": gs.get("n_hits", 0),
                    "best_guide_length": gs.get("best_length", 0),
                })
        logger.info("Wrote summary to %s", tsv_path)

    # ------------------------------------------------------------------
    # All-in-one
    # ------------------------------------------------------------------

    def run(
        self,
        protein_fasta: str,
        formatter_dir: str,
        output_dir: str,
        cpus: int = 8,
        min_tail_head: int = 1,
    ) -> List[Dict]:
        """Full pipeline: hmmsearch → filter → export.

        Args:
            protein_fasta: path to all_proteins.faa from system clustering.
            formatter_dir: directory with sample subdirs containing guide JSONs.
            output_dir: directory for output files.
            cpus: threads for hmmsearch.
            min_tail_head: minimum tail-head reads to pass.

        Returns:
            List of matching record dicts.
        """
        is110_ids, trans_domains = self.run_hmmsearch(protein_fasta, output_dir, cpus=cpus)
        records = self.filter_records(
            formatter_dir, is110_ids,
            min_tail_head=min_tail_head,
            trans_domains=trans_domains,
        )
        self.export_results(records, output_dir)
        return records
