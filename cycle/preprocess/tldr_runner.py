"""Run tldr per organism group to detect transposon insertions."""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import pandas as pd

from ..utils import slugify as _slugify
from .config import DEFAULT_ALIGNMENT_DIR, DEFAULT_TLDR_OUTPUT_DIR

logger = logging.getLogger(__name__)


class TldrRunner:
    """Run tldr once per organism group.

    tldr requires all BAMs in a single invocation to share the same
    reference genome, so samples are grouped by organism before calling.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_TLDR_OUTPUT_DIR,
        alignment_dir: str = DEFAULT_ALIGNMENT_DIR,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.alignment_dir = Path(alignment_dir)

        if not shutil.which("tldr"):
            raise RuntimeError("tldr not found in PATH")

    def run_organism_group(
        self,
        bams: list[Path],
        ref_fasta: Path,
        is_ref: Path,
        organism: str,
        procs: int = 8,
    ) -> Optional[Path]:
        """Run tldr on a group of BAMs sharing the same reference genome.

        Args:
            bams: Sorted BAM files (all aligned to the same reference).
            ref_fasta: Reference genome FASTA used for alignment.
            is_ref: tldr-formatted IS element reference FASTA.
            organism: Organism name (used for output directory naming).
            procs: Number of processes for tldr.

        Returns:
            Path to the tldr output table, or None on failure.
        """
        slug = _slugify(organism)
        org_dir = self.output_dir / slug
        org_dir.mkdir(parents=True, exist_ok=True)

        output_prefix = org_dir / slug
        table_path = Path(f"{output_prefix}.table.txt")

        if table_path.exists():
            logger.info(f"tldr output exists for {organism}: {table_path}")
            return table_path

        cmd = [
            "tldr",
            "-b", ",".join(str(b) for b in bams),
            "-e", str(is_ref),
            "-r", str(ref_fasta),
            "-p", str(procs),
            "-o", str(output_prefix),
        ]

        logger.info(
            f"Running tldr for {organism} ({len(bams)} BAMs): {slug}"
        )
        logger.debug(f"  cmd: {' '.join(cmd)}")

        try:
            ret = subprocess.run(
                cmd, capture_output=True, text=True, timeout=7200,
                cwd=str(org_dir),
            )
        except subprocess.TimeoutExpired:
            logger.error(f"tldr timed out for {organism}")
            return None

        if ret.returncode != 0:
            logger.error(
                f"tldr failed for {organism} (rc={ret.returncode}): "
                f"{ret.stderr.strip()[:500]}"
            )
            return None

        if table_path.exists():
            logger.info(f"  -> {table_path}")
            return table_path

        logger.warning(f"tldr finished but no table found at {table_path}")
        return None

    def run_batch(
        self,
        metadata: pd.DataFrame,
        ref_map: dict[str, dict],
        is_ref: Path,
        procs: int = 8,
        organism_col: str = "organism",
        accession_col: str = "srr_accession",
    ) -> dict[str, Optional[Path]]:
        """Group samples by organism and run tldr for each group.

        Args:
            metadata: DataFrame with organism and accession columns.
            ref_map: Organism -> {accession, fasta, fai} from resolve step.
            is_ref: Path to IS element reference FASTA.
            procs: Number of processes for tldr.
            organism_col: Column name for organism.
            accession_col: Column name for SRR accession.

        Returns:
            Dict mapping organism -> tldr table path (or None).
        """
        results: dict[str, Optional[Path]] = {}
        groups = metadata.groupby(organism_col)

        logger.info(f"Running tldr for {len(groups)} organism groups")

        for organism, group_df in groups:
            ref_info = ref_map.get(organism)
            if not ref_info:
                logger.warning(f"No reference for {organism}, skipping tldr")
                results[organism] = None
                continue

            ref_fasta = Path(ref_info["fasta"])

            # Collect BAMs for this organism
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

            table = self.run_organism_group(
                bams=bams,
                ref_fasta=ref_fasta,
                is_ref=is_ref,
                organism=organism,
                procs=procs,
            )
            results[organism] = table

        ok = sum(1 for v in results.values() if v)
        logger.info(f"tldr complete: {ok}/{len(results)} organism groups succeeded")
        return results
