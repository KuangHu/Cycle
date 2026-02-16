"""Align nanopore reads to reference genomes with minimap2 → sorted BAM."""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from .config import (
    DEFAULT_ALIGNMENT_DIR,
    DEFAULT_MINIMAP2_PRESET,
    DEFAULT_SORT_MEMORY,
    DEFAULT_THREADS,
)

logger = logging.getLogger(__name__)


class Aligner:
    """Minimap2 index + align → sorted BAM + .bai.

    One ``.mmi`` index is built per reference genome (shared across samples
    of the same organism).  Alignment pipes minimap2 directly into
    ``samtools sort`` with no intermediate SAM on disk.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_ALIGNMENT_DIR,
        preset: str = DEFAULT_MINIMAP2_PRESET,
        threads: int = DEFAULT_THREADS,
        sort_memory: str = DEFAULT_SORT_MEMORY,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.preset = preset
        self.threads = threads
        self.sort_memory = sort_memory

        for tool in ("minimap2", "samtools"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def index(self, fasta: Path) -> Optional[Path]:
        """Build a minimap2 ``.mmi`` index for a reference FASTA.

        Skips if the index already exists. Returns the path to the ``.mmi``
        file, or None on failure.
        """
        fasta = Path(fasta)
        mmi = fasta.with_suffix(".mmi")

        if mmi.exists():
            logger.info(f"Index exists: {mmi}")
            return mmi

        cmd = [
            "minimap2",
            "-x", self.preset,
            "-t", str(self.threads),
            "-d", str(mmi),
            str(fasta),
        ]
        logger.info(f"Indexing {fasta.name}")
        ret = subprocess.run(cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            logger.error(f"minimap2 index failed: {ret.stderr.strip()}")
            return None

        logger.info(f"  → {mmi}")
        return mmi

    def align(
        self,
        fastq: Path,
        reference: Path,
        sample_id: Optional[str] = None,
    ) -> Optional[Path]:
        """Align reads to a reference and produce a sorted BAM + index.

        Args:
            fastq: Path to input FASTQ (or .fastq.gz).
            reference: Path to reference FASTA or ``.mmi`` index.
            sample_id: Used for output filename. Defaults to FASTQ stem.

        Returns:
            Path to the sorted BAM, or None on failure.
        """
        fastq = Path(fastq)
        reference = Path(reference)
        if sample_id is None:
            sample_id = fastq.name.split(".")[0]

        bam = self.output_dir / f"{sample_id}.sorted.bam"

        if bam.exists() and Path(str(bam) + ".bai").exists():
            logger.info(f"BAM exists: {bam}")
            return bam

        # Prefer .mmi if available
        mmi = reference.with_suffix(".mmi")
        ref = mmi if mmi.exists() else reference

        # minimap2 -a | samtools sort → BAM (pipe-based, no intermediate SAM)
        # --MD flag is required for Sniffles2 to parse alignments correctly
        mm2_cmd = [
            "minimap2",
            "-a",
            "--MD",
            "-x", self.preset,
            "-t", str(self.threads),
            str(ref),
            str(fastq),
        ]
        sort_cmd = [
            "samtools", "sort",
            "-@", str(self.threads),
            "-m", self.sort_memory,
            "-o", str(bam),
            "-",
        ]

        logger.info(f"Aligning {fastq.name} → {bam.name}")
        try:
            mm2_proc = subprocess.Popen(
                mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            sort_proc = subprocess.Popen(
                sort_cmd,
                stdin=mm2_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            # Allow mm2 to receive SIGPIPE if sort exits early
            mm2_proc.stdout.close()

            sort_out, sort_err = sort_proc.communicate()
            mm2_err = mm2_proc.stderr.read()
            mm2_proc.stderr.close()
            mm2_proc.wait()

            if sort_proc.returncode != 0:
                logger.error(
                    f"samtools sort failed for {sample_id}: "
                    f"{sort_err.decode().strip()}"
                )
                return None

        except Exception as exc:
            logger.error(f"Alignment pipeline failed for {sample_id}: {exc}")
            return None

        # Index the BAM
        ret = subprocess.run(
            ["samtools", "index", str(bam)],
            capture_output=True, text=True,
        )
        if ret.returncode != 0:
            logger.error(f"samtools index failed: {ret.stderr.strip()}")
            return None

        logger.info(f"  → {bam} + .bai")
        return bam

    def align_batch(
        self,
        sample_map: list[dict],
    ) -> list[dict]:
        """Align a batch of samples.

        Args:
            sample_map: List of dicts with keys
                {sample_id, fastq, reference_fasta}.

        Returns:
            List of dicts with added ``bam`` key (Path or None).
        """
        results = []
        for i, entry in enumerate(sample_map, 1):
            sid = entry["sample_id"]
            fastq = Path(entry["fastq"])
            ref = Path(entry["reference_fasta"])

            logger.info(
                f"[{i}/{len(sample_map)}] {sid}"
            )

            # Ensure index exists
            self.index(ref)

            bam = self.align(fastq, ref, sample_id=sid)
            results.append({**entry, "bam": bam})

        ok = sum(1 for r in results if r["bam"])
        logger.info(f"Alignment complete: {ok}/{len(results)} succeeded")
        return results
