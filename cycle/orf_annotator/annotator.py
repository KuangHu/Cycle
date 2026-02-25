"""ORF annotation for IS element sequences using Prodigal."""

from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class ORFAnnotator:
    """Predict ORFs on IS element sequences with Prodigal (meta mode)."""

    def __init__(self) -> None:
        if not shutil.which("prodigal"):
            raise RuntimeError("prodigal not found in PATH")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def annotate_sample(self, json_path: str | Path) -> Optional[Path]:
        """Annotate IS records in a formatter JSON with ORF predictions.

        Parameters
        ----------
        json_path : path to ``*_is_records.json`` from the IS formatter.

        Returns
        -------
        Path to the annotated JSON (``*_is_records_annotated.json``),
        or None on failure.
        """
        json_path = Path(json_path)
        if not json_path.exists():
            logger.warning(f"File not found: {json_path}")
            return None

        with open(json_path) as fh:
            records = json.load(fh)

        if not records:
            logger.warning(f"No records in {json_path.name}")
            return None

        # Collect sequences keyed by index for batch prodigal call
        seqs: dict[int, str] = {}
        for idx, rec in enumerate(records):
            seq = (rec.get("is_element") or {}).get("sequence", "")
            if seq:
                seqs[idx] = seq

        if not seqs:
            logger.info(f"No IS sequences to annotate in {json_path.name}")
            for rec in records:
                rec["orf_annotation"] = None
            return self._write_output(json_path, records)

        # Run prodigal once on all sequences
        gff_text, protein_text = self._run_prodigal(seqs)
        if gff_text is None:
            logger.error(f"Prodigal failed for {json_path.name}")
            return None

        # Parse results per sequence
        orfs_by_idx = self._parse_gff(gff_text)
        proteins_by_key = self._parse_proteins(protein_text)

        # Attach annotations
        for idx, rec in enumerate(records):
            if idx not in seqs:
                rec["orf_annotation"] = None
                continue

            seq_len = len(seqs[idx])
            orfs = orfs_by_idx.get(idx, [])

            # Attach protein sequences from FASTA
            for orf in orfs:
                key = (idx, orf["start"], orf["end"], orf["strand"])
                orf["protein_sequence"] = proteins_by_key.get(key, "")

            noncoding = self._compute_noncoding(orfs, seq_len)
            coding_nt = sum(o["length_nt"] for o in orfs)

            rec["orf_annotation"] = {
                "num_orfs": len(orfs),
                "coding_fraction": round(coding_nt / seq_len, 4) if seq_len else 0,
                "orfs": orfs,
                "noncoding_regions": noncoding,
            }

        sample_id = json_path.stem.replace("_is_records", "")
        n_annotated = sum(1 for r in records if r.get("orf_annotation"))
        total_orfs = sum(
            r["orf_annotation"]["num_orfs"]
            for r in records
            if r.get("orf_annotation")
        )
        logger.info(
            f"{sample_id}: annotated {n_annotated}/{len(records)} records, "
            f"{total_orfs} ORFs total"
        )

        return self._write_output(json_path, records)

    def annotate_batch(
        self, formatter_dir: str | Path, parallel: int = 1
    ) -> dict[str, Optional[Path]]:
        """Annotate all ``*_is_records.json`` files in *formatter_dir*.

        Returns dict mapping sample JSON filename to output path (or None).
        """
        formatter_dir = Path(formatter_dir)
        json_files = sorted(formatter_dir.glob("*/*_is_records.json"))

        if not json_files:
            logger.warning(f"No *_is_records.json files in {formatter_dir}")
            return {}

        logger.info(
            f"Found {len(json_files)} sample JSON files in {formatter_dir}"
        )

        results: dict[str, Optional[Path]] = {}

        if parallel <= 1:
            for jf in json_files:
                try:
                    results[jf.name] = self.annotate_sample(jf)
                except Exception as exc:
                    logger.error(f"Failed for {jf.name}: {exc}")
                    results[jf.name] = None
        else:
            from concurrent.futures import ProcessPoolExecutor, as_completed

            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for jf in json_files:
                    fut = pool.submit(_annotate_worker, jf)
                    futures[fut] = jf.name

                for fut in as_completed(futures):
                    name = futures[fut]
                    try:
                        results[name] = fut.result()
                    except Exception as exc:
                        logger.error(f"Worker failed for {name}: {exc}")
                        results[name] = None

        succeeded = sum(1 for v in results.values() if v is not None)
        logger.info(f"Batch done: {succeeded}/{len(results)} samples annotated")
        return results

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _run_prodigal(
        self, seqs: dict[int, str]
    ) -> tuple[Optional[str], Optional[str]]:
        """Run prodigal -p meta on batched sequences.

        Returns (gff_text, protein_fasta_text) or (None, None) on failure.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "input.fasta"
            gff_path = Path(tmpdir) / "output.gff"
            prot_path = Path(tmpdir) / "output.faa"

            # Write multi-FASTA; use index as sequence ID
            with open(fasta_path, "w") as fh:
                for idx, seq in seqs.items():
                    fh.write(f">seq_{idx}\n{seq}\n")

            cmd = [
                "prodigal",
                "-i", str(fasta_path),
                "-o", str(gff_path),
                "-a", str(prot_path),
                "-f", "gff",
                "-p", "meta",
                "-q",
            ]

            ret = subprocess.run(cmd, capture_output=True, text=True)
            if ret.returncode != 0:
                logger.error(f"prodigal failed: {ret.stderr.strip()[:500]}")
                return None, None

            gff_text = gff_path.read_text() if gff_path.exists() else ""
            prot_text = prot_path.read_text() if prot_path.exists() else ""
            return gff_text, prot_text

    def _parse_gff(self, gff_text: str) -> dict[int, list[dict]]:
        """Parse prodigal GFF output into ORF dicts keyed by sequence index."""
        orfs_by_idx: dict[int, list[dict]] = {}

        for line in gff_text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue

            # Sequence ID is "seq_{idx}"
            m = re.match(r"seq_(\d+)", parts[0])
            if not m:
                continue
            idx = int(m.group(1))

            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            score = parts[5]

            # Parse confidence from attributes (conf=X.XX)
            attrs = parts[8]
            conf_match = re.search(r"conf=([0-9.]+)", attrs)
            confidence = float(conf_match.group(1)) if conf_match else None

            orf = {
                "start": start,
                "end": end,
                "strand": strand,
                "length_nt": end - start + 1,
                "protein_sequence": "",  # filled later
                "confidence": confidence,
            }

            orfs_by_idx.setdefault(idx, []).append(orf)

        # Sort ORFs by start position within each sequence
        for idx in orfs_by_idx:
            orfs_by_idx[idx].sort(key=lambda o: o["start"])

        return orfs_by_idx

    def _parse_proteins(self, protein_text: str) -> dict[tuple, str]:
        """Parse prodigal protein FASTA.

        Returns dict keyed by (seq_idx, start, end, strand) -> protein seq.
        """
        proteins: dict[tuple, str] = {}
        current_key: Optional[tuple] = None
        current_seq: list[str] = []

        for line in protein_text.splitlines():
            if line.startswith(">"):
                # Save previous
                if current_key is not None:
                    proteins[current_key] = "".join(current_seq)

                # Parse header: >seq_0_1 # 1 # 300 # 1 # ID=...
                parts = line[1:].split(" # ")
                if len(parts) >= 4:
                    header_id = parts[0].strip()
                    m = re.match(r"seq_(\d+)_\d+", header_id)
                    if m:
                        idx = int(m.group(1))
                        start = int(parts[1].strip())
                        end = int(parts[2].strip())
                        strand = "+" if parts[3].strip() == "1" else "-"
                        current_key = (idx, start, end, strand)
                    else:
                        current_key = None
                else:
                    current_key = None
                current_seq = []
            else:
                current_seq.append(line.strip())

        # Save last entry
        if current_key is not None:
            proteins[current_key] = "".join(current_seq)

        return proteins

    @staticmethod
    def _compute_noncoding(orfs: list[dict], seq_len: int) -> list[dict]:
        """Compute noncoding regions (gaps between CDS features)."""
        if not orfs:
            return [{"start": 1, "end": seq_len, "length": seq_len, "type": "full"}]

        regions: list[dict] = []
        first_start = orfs[0]["start"]
        last_end = orfs[-1]["end"]

        # 5-prime noncoding
        if first_start > 1:
            length = first_start - 1
            regions.append({
                "start": 1,
                "end": first_start - 1,
                "length": length,
                "type": "5_prime",
            })

        # Intergenic regions
        for i in range(len(orfs) - 1):
            gap_start = orfs[i]["end"] + 1
            gap_end = orfs[i + 1]["start"] - 1
            if gap_end >= gap_start:
                regions.append({
                    "start": gap_start,
                    "end": gap_end,
                    "length": gap_end - gap_start + 1,
                    "type": "intergenic",
                })

        # 3-prime noncoding
        if last_end < seq_len:
            length = seq_len - last_end
            regions.append({
                "start": last_end + 1,
                "end": seq_len,
                "length": length,
                "type": "3_prime",
            })

        return regions

    @staticmethod
    def _write_output(json_path: Path, records: list[dict]) -> Path:
        """Write annotated records to ``*_is_records_annotated.json``."""
        out_path = json_path.with_name(
            json_path.name.replace("_is_records.json", "_is_records_annotated.json")
        )
        with open(out_path, "w") as fh:
            json.dump(records, fh, indent=2)
        logger.info(f"Wrote {len(records)} annotated records to {out_path}")
        return out_path


# ------------------------------------------------------------------
# Worker for parallel execution
# ------------------------------------------------------------------

def _annotate_worker(json_path: Path) -> Optional[Path]:
    """Standalone worker for ProcessPoolExecutor."""
    annotator = ORFAnnotator()
    return annotator.annotate_sample(json_path)
