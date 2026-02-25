"""
IS110 Circular Element Filter

Identifies IS elements that:
  1. Have tail-head junction reads (circle_evidence.n_tail_head_reads > 0)
  2. Contain IS110-family transposase (DEDD or Tnp20 HMM hit)

Workflow:
  1. Run hmmsearch with DEDD.hmm and/or Tnp20.hmm against a protein FASTA
  2. Parse hits → set of transposon IDs containing IS110 protein
  3. Scan *_is_records_guide.json files, emit records matching both criteria
"""

import csv
import json
import logging
import os
import subprocess
import tempfile
from glob import glob
from typing import Dict, List, Optional, Set

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
    ):
        self.dedd_hmm = dedd_hmm
        self.tnp20_hmm = tnp20_hmm
        self.evalue = evalue

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
    ) -> Set[str]:
        """Run both DEDD and Tnp20 hmmsearch, return set of transposon IDs with IS110 hits."""
        os.makedirs(output_dir, exist_ok=True)
        is110_ids = set()

        for name, hmm in [("DEDD", self.dedd_hmm), ("Tnp20", self.tnp20_hmm)]:
            tbl = os.path.join(output_dir, f"{name}_hits.tbl")
            self.hmmsearch(protein_fasta, tbl, hmm, cpus=cpus)
            ids = self.parse_tblout(tbl)
            logger.info("%s: %d protein hits → %d unique transposons", name, self._last_n_hits, len(ids))
            is110_ids.update(ids)

        logger.info("Combined: %d unique IS110 transposons", len(is110_ids))
        return is110_ids

    # ------------------------------------------------------------------
    # Step 2: parse hmmsearch output
    # ------------------------------------------------------------------

    def parse_tblout(self, tbl_path: str) -> Set[str]:
        """Parse hmmsearch tblout, extract transposon IDs from protein headers.

        Header format: {transposon_id}__{start}_{end}_{strand}
        """
        ids = set()
        n_hits = 0
        with open(tbl_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                n_hits += 1
                protein_id = line.split()[0]
                # Split on __ to get transposon ID
                transposon_id = protein_id.rsplit("__", 1)[0]
                ids.add(transposon_id)
        self._last_n_hits = n_hits
        return ids

    # ------------------------------------------------------------------
    # Step 3: filter records
    # ------------------------------------------------------------------

    def filter_records(
        self,
        formatter_dir: str,
        is110_ids: Set[str],
        min_tail_head: int = 1,
    ) -> List[Dict]:
        """Scan guide JSONs, return records with IS110 protein AND tail-head reads.

        Args:
            formatter_dir: directory with sample subdirs containing *_is_records_guide.json.
            is110_ids: set of transposon IDs that contain IS110 protein.
            min_tail_head: minimum n_tail_head_reads to pass filter.

        Returns:
            List of matching record dicts.
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
                    matched.append(rec)

        logger.info(
            "Scanned %d records: %d IS110, %d circular (TH>=%d), %d both",
            total_scanned, total_is110, total_circular, min_tail_head, len(matched),
        )
        return matched

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
        is110_ids = self.run_hmmsearch(protein_fasta, output_dir, cpus=cpus)
        records = self.filter_records(formatter_dir, is110_ids, min_tail_head=min_tail_head)
        self.export_results(records, output_dir)
        return records
