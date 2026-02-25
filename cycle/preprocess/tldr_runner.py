"""Run tldr per sample to detect transposon insertions."""

import logging
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pandas as pd

from .config import DEFAULT_ALIGNMENT_DIR, DEFAULT_TLDR_OUTPUT_DIR

logger = logging.getLogger(__name__)


class TldrRunner:
    """Run tldr once per sample."""

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

    def run_sample(
        self,
        bam: Path,
        ref_fasta: Path,
        is_ref: Path,
        sample_id: str,
        procs: int = 8,
    ) -> Optional[Path]:
        """Run tldr on a single BAM for one sample.

        Args:
            bam: Sorted BAM file.
            ref_fasta: Reference genome FASTA used for alignment.
            is_ref: tldr-formatted IS element reference FASTA.
            sample_id: Sample accession (used for output directory naming).
            procs: Number of processes for tldr.

        Returns:
            Path to the tldr output table, or None on failure.
        """
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        output_prefix = sample_dir / sample_id
        table_path = Path(f"{output_prefix}.table.txt")

        if table_path.exists():
            logger.info(f"tldr output exists for {sample_id}: {table_path}")
            return table_path

        cmd = [
            "tldr",
            "-b", str(bam),
            "-e", str(is_ref),
            "-r", str(ref_fasta),
            "-p", str(procs),
            "-o", str(output_prefix),
        ]

        logger.info(f"Running tldr for {sample_id}: {bam.name}")
        logger.debug(f"  cmd: {' '.join(cmd)}")

        try:
            ret = subprocess.run(
                cmd, capture_output=True, text=True, timeout=7200,
                cwd=str(sample_dir),
            )
        except subprocess.TimeoutExpired:
            logger.error(f"tldr timed out for {sample_id}")
            return None

        if ret.returncode != 0:
            logger.error(
                f"tldr failed for {sample_id} (rc={ret.returncode}): "
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
        ref_map: dict[str, Path],
        is_ref: Path,
        procs: int = 8,
        parallel: int = 1,
        accession_col: str = "srr_accession",
    ) -> dict[str, Optional[Path]]:
        """Run tldr for each sample.

        Args:
            metadata: DataFrame with accession column.
            ref_map: sample_id -> ref_fasta path.
            is_ref: Path to IS element reference FASTA.
            procs: Number of processes per tldr invocation.
            parallel: Number of samples to run in parallel.
            accession_col: Column name for SRR accession.

        Returns:
            Dict mapping sample_id -> tldr table path (or None).
        """
        results: dict[str, Optional[Path]] = {}

        # Build task list
        tasks: list[tuple[str, Path, Path]] = []
        for _, row in metadata.iterrows():
            sid = row[accession_col]
            ref_fasta = ref_map.get(sid)
            if not ref_fasta:
                logger.warning(f"No reference for {sid}, skipping tldr")
                results[sid] = None
                continue

            bam = self.alignment_dir / f"{sid}.sorted.bam"
            if not bam.exists():
                logger.warning(f"BAM not found for {sid}: {bam}")
                results[sid] = None
                continue

            tasks.append((sid, bam, ref_fasta))

        logger.info(f"Running tldr for {len(tasks)} samples")

        if parallel <= 1:
            for sid, bam, ref_fasta in tasks:
                table = self.run_sample(
                    bam=bam, ref_fasta=ref_fasta, is_ref=is_ref,
                    sample_id=sid, procs=procs,
                )
                results[sid] = table
        else:
            logger.info(f"Running {len(tasks)} samples with {parallel} in parallel, {procs} procs each")
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for sid, bam, ref_fasta in tasks:
                    fut = pool.submit(
                        _run_tldr_worker,
                        bam=bam,
                        ref_fasta=ref_fasta,
                        is_ref=is_ref,
                        sample_id=sid,
                        procs=procs,
                        output_dir=self.output_dir,
                    )
                    futures[fut] = sid

                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        table = fut.result()
                        results[sid] = table
                        if table:
                            logger.info(f"  -> {table}")
                    except Exception as e:
                        logger.error(f"tldr worker failed for {sid}: {e}")
                        results[sid] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(f"tldr complete: {ok}/{len(results)} samples succeeded")
        return results


def _run_tldr_worker(
    bam: Path,
    ref_fasta: Path,
    is_ref: Path,
    sample_id: str,
    procs: int,
    output_dir: Path,
) -> Optional[Path]:
    """Standalone worker function for parallel tldr execution."""
    sample_dir = output_dir / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    output_prefix = sample_dir / sample_id
    table_path = Path(f"{output_prefix}.table.txt")

    if table_path.exists():
        return table_path

    cmd = [
        "tldr",
        "-b", str(bam),
        "-e", str(is_ref),
        "-r", str(ref_fasta),
        "-p", str(procs),
        "-o", str(output_prefix),
    ]

    try:
        ret = subprocess.run(
            cmd, capture_output=True, text=True, timeout=7200,
            cwd=str(sample_dir),
        )
    except subprocess.TimeoutExpired:
        logging.getLogger(__name__).error(f"tldr timed out for {sample_id}")
        return None

    if ret.returncode != 0:
        logging.getLogger(__name__).error(
            f"tldr failed for {sample_id} (rc={ret.returncode}): "
            f"{ret.stderr.strip()[:500]}"
        )
        return None

    if table_path.exists():
        return table_path

    return None
