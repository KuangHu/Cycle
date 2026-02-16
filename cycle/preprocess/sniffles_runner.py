"""Run Sniffles2 per organism to detect insertions as IS candidates."""

import logging
import shutil
import subprocess
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pandas as pd

from ..utils import slugify as _slugify
from .config import DEFAULT_ALIGNMENT_DIR

logger = logging.getLogger(__name__)


class SnifflesRunner:
    """Run Sniffles2 once per organism group to detect insertions.

    Sniffles2 is a structural variant caller that can detect insertions
    without requiring a transposon reference library. We use it as a
    faster alternative to tldr for finding IS insertion candidates.
    """

    def __init__(
        self,
        output_dir: str,
        alignment_dir: str = DEFAULT_ALIGNMENT_DIR,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.alignment_dir = Path(alignment_dir)

        if not shutil.which("sniffles"):
            raise RuntimeError("sniffles not found in PATH")

    def run_organism_group(
        self,
        bams: list[Path],
        ref_fasta: Path,
        organism: str,
        min_size: int = 500,  # Lower threshold to catch smaller IS elements
        max_size: int = 20000,  # Upper limit to include larger IS elements and composite transposons
        min_support: int = 1,  # Detect rare transposition events with 1 read
        disable_qc: bool = True,  # Disable QC filters for max sensitivity
    ) -> Optional[Path]:
        """Run Sniffles2 on a group of BAMs for one organism.

        Args:
            bams: Sorted BAM files (all aligned to the same reference).
            ref_fasta: Reference genome FASTA used for alignment.
            organism: Organism name (used for output directory naming).
            min_size: Minimum insertion size in bp (IS elements typically 500-3000bp).
            max_size: Maximum insertion size in bp (includes composite transposons).
            min_support: Minimum number of supporting reads.
            disable_qc: Disable Sniffles2 quality control filters for max sensitivity.

        Returns:
            Path to the Sniffles2-derived insertion table, or None on failure.
        """
        slug = _slugify(organism)
        org_dir = self.output_dir / slug
        org_dir.mkdir(parents=True, exist_ok=True)

        output_prefix = org_dir / slug
        table_path = Path(f"{output_prefix}.table.txt")
        vcf_path = Path(f"{output_prefix}.vcf")

        if table_path.exists():
            logger.info(f"Sniffles2 output exists for {organism}: {table_path}")
            return table_path

        # Run Sniffles2 with sensitive settings for detecting rare transpositions
        # Note: --mosaic flag is NOT used because it's designed for somatic/cancer
        # variants and filters out bacterial transposon insertions (tested on
        # A. xylosoxidans: 2 insertions without --mosaic, 0 with --mosaic).
        # --no-qc disables quality filters for maximum sensitivity (found 235
        # insertions vs 2 with default QC on A. xylosoxidans test).
        cmd = [
            "sniffles",
            "--input", ",".join(str(b) for b in bams),
            "--vcf", str(vcf_path),
            "--reference", str(ref_fasta),
            "--threads", "4",
            "--minsvlen", str(min_size),
            "--minsupport", str(min_support),
        ]
        if disable_qc:
            cmd.append("--no-qc")

        logger.info(
            f"Running Sniffles2 for {organism} ({len(bams)} BAMs): {slug}"
        )
        logger.debug(f"  cmd: {' '.join(cmd)}")

        try:
            ret = subprocess.run(
                cmd, capture_output=True, text=True, timeout=1800,
                cwd=str(org_dir),
            )
        except subprocess.TimeoutExpired:
            logger.error(f"Sniffles2 timed out for {organism}")
            return None

        if ret.returncode != 0:
            logger.error(
                f"Sniffles2 failed for {organism} (rc={ret.returncode}): "
                f"{ret.stderr.strip()[:500]}"
            )
            return None

        if not vcf_path.exists():
            logger.warning(f"Sniffles2 finished but no VCF found at {vcf_path}")
            return None

        # Parse VCF and create insertion table
        try:
            insertions = self._parse_vcf(vcf_path, min_size, max_size, organism)
            if insertions:
                self._write_table(insertions, table_path)
                logger.info(f"  -> {table_path} ({len(insertions)} insertions)")
                return table_path
            else:
                logger.warning(f"No insertions found in {vcf_path}")
                # Write empty table
                pd.DataFrame(columns=["UUID", "Chrom", "Start", "End", "Sequence", "Support"]).to_csv(
                    table_path, sep="\t", index=False
                )
                return table_path
        except Exception as e:
            logger.error(f"Failed to parse VCF for {organism}: {e}")
            return None

    def _parse_vcf(
        self, vcf_path: Path, min_size: int, max_size: int, organism: str
    ) -> list[dict]:
        """Parse Sniffles2 VCF and extract insertions."""
        insertions = []

        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                # Parse INFO field for SVTYPE and SVLEN
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, val = item.split("=", 1)
                        info_dict[key] = val

                # Filter for insertions, duplications, and breakends (all relevant for IS detection)
                svtype = info_dict.get("SVTYPE", "")
                if svtype not in ("INS", "DUP", "BND"):
                    continue

                # Get variant length
                svlen = info_dict.get("SVLEN")
                if svlen:
                    svlen = abs(int(svlen))
                    if svlen < min_size or svlen > max_size:
                        continue
                elif svtype == "BND":
                    # Breakends don't have length, treat as min_size for now
                    svlen = min_size
                else:
                    # Try to infer from ALT
                    if alt.startswith("<"):
                        # Symbolic alt without length info, skip if no SVLEN
                        continue
                    svlen = len(alt) - len(ref)
                    if svlen < min_size or svlen > max_size:
                        continue

                # Get support count
                support = int(info_dict.get("SUPPORT", "0"))

                # Extract inserted sequence from ALT if available
                if alt.startswith("<"):
                    # Symbolic alt, no sequence
                    sequence = "N" * svlen
                else:
                    # Sequence-resolved insertion
                    sequence = alt[len(ref):]

                insertions.append({
                    "UUID": str(uuid.uuid4()),
                    "Chrom": chrom,
                    "Start": pos,
                    "End": pos + svlen,
                    "Sequence": sequence,
                    "Support": support,
                    "Family": f"Sniffles2_{svtype}",  # Tag with variant type
                    "Subfamily": f"{organism}_{chrom}_{pos}",
                })

        return insertions

    def _write_table(self, insertions: list[dict], output_path: Path):
        """Write insertion table in tldr-compatible format."""
        df = pd.DataFrame(insertions)
        # Rename Sequence to Consensus to match tldr format
        df = df.rename(columns={"Sequence": "Consensus"})
        # Reorder columns to match tldr format expectations
        cols = ["UUID", "Chrom", "Start", "End", "Family", "Subfamily", "Consensus", "Support"]
        df = df[cols]
        df.to_csv(output_path, sep="\t", index=False)

    def run_batch(
        self,
        metadata: pd.DataFrame,
        ref_map: dict[str, dict],
        parallel: int = 1,
        min_size: int = 500,
        max_size: int = 20000,
        disable_qc: bool = True,
        organism_col: str = "organism",
        accession_col: str = "srr_accession",
    ) -> dict[str, Optional[Path]]:
        """Group samples by organism and run Sniffles2 for each group.

        Args:
            metadata: DataFrame with organism and accession columns.
            ref_map: Organism -> {accession, fasta, fai} from resolve step.
            parallel: Number of organisms to run in parallel.
            min_size: Minimum insertion size in bp.
            max_size: Maximum insertion size in bp.
            disable_qc: Disable Sniffles2 quality control filters.
            organism_col: Column name for organism.
            accession_col: Column name for SRR accession.

        Returns:
            Dict mapping organism -> insertion table path (or None).
        """
        results: dict[str, Optional[Path]] = {}
        groups = metadata.groupby(organism_col)

        logger.info(f"Running Sniffles2 for {len(groups)} organism groups")

        # Build task list
        tasks: list[tuple[str, list[Path], Path]] = []
        for organism, group_df in groups:
            ref_info = ref_map.get(organism)
            if not ref_info:
                logger.warning(f"No reference for {organism}, skipping Sniffles2")
                results[organism] = None
                continue

            ref_fasta = Path(ref_info["fasta"])

            bams = []
            for _, row in group_df.iterrows():
                sid = row[accession_col]
                bam = self.alignment_dir / f"{sid}.sorted.bam"
                if bam.exists():
                    bams.append(bam)
                else:
                    logger.warning(f"BAM not found for {sid}: {bam}")

            if not bams:
                logger.warning(f"No BAMs found for {organism}, skipping")
                results[organism] = None
                continue

            tasks.append((organism, bams, ref_fasta))

        if parallel <= 1:
            # Sequential
            for organism, bams, ref_fasta in tasks:
                table = self.run_organism_group(
                    bams=bams, ref_fasta=ref_fasta,
                    organism=organism, min_size=min_size, max_size=max_size,
                    disable_qc=disable_qc,
                )
                results[organism] = table
        else:
            # Parallel execution
            logger.info(f"Running {len(tasks)} organisms with {parallel} in parallel")
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for organism, bams, ref_fasta in tasks:
                    fut = pool.submit(
                        _run_sniffles_worker,
                        bams=bams,
                        ref_fasta=ref_fasta,
                        organism=organism,
                        min_size=min_size,
                        max_size=max_size,
                        disable_qc=disable_qc,
                        output_dir=self.output_dir,
                        alignment_dir=self.alignment_dir,
                    )
                    futures[fut] = organism

                for fut in as_completed(futures):
                    organism = futures[fut]
                    try:
                        table = fut.result()
                        results[organism] = table
                        if table:
                            logger.info(f"  -> {table}")
                    except Exception as e:
                        logger.error(f"Sniffles2 worker failed for {organism}: {e}")
                        results[organism] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(f"Sniffles2 complete: {ok}/{len(results)} organism groups succeeded")
        return results


def _run_sniffles_worker(
    bams: list[Path],
    ref_fasta: Path,
    organism: str,
    min_size: int,
    max_size: int,
    disable_qc: bool,
    output_dir: Path,
    alignment_dir: Path,
) -> Optional[Path]:
    """Standalone worker function for parallel Sniffles2 execution."""
    runner = SnifflesRunner(
        output_dir=str(output_dir),
        alignment_dir=str(alignment_dir),
    )
    return runner.run_organism_group(
        bams=bams,
        ref_fasta=ref_fasta,
        organism=organism,
        min_size=min_size,
        max_size=max_size,
        disable_qc=disable_qc,
    )
