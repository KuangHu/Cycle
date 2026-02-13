"""Download SRA runs using kingfisher (multi-source with automatic fallback)."""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import Optional

import pandas as pd

from .config import (
    DEFAULT_DOWNLOAD_DIR,
    DEFAULT_DOWNLOAD_METHODS,
    DEFAULT_OUTPUT_FORMAT,
    DEFAULT_DOWNLOAD_THREADS,
    DEFAULT_EXTRACTION_THREADS,
)

logger = logging.getLogger(__name__)


class SRADownloader:
    """Batch download SRA runs to FASTQ using kingfisher."""

    def __init__(
        self,
        output_dir: str = DEFAULT_DOWNLOAD_DIR,
        methods: Optional[list[str]] = None,
        output_format: str = DEFAULT_OUTPUT_FORMAT,
        download_threads: int = DEFAULT_DOWNLOAD_THREADS,
        extraction_threads: int = DEFAULT_EXTRACTION_THREADS,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.methods = methods or DEFAULT_DOWNLOAD_METHODS
        self.output_format = output_format
        self.download_threads = download_threads
        self.extraction_threads = extraction_threads

        if not shutil.which("kingfisher"):
            raise RuntimeError("kingfisher not found in PATH. Install kingfisher.")

    def _build_cmd(self, accessions: list[str]) -> list[str]:
        """Build the kingfisher get command."""
        cmd = [
            "kingfisher", "get",
            "-r", *accessions,
            "-m", *self.methods,
            "-f", self.output_format,
            "--output-directory", str(self.output_dir),
            "--download-threads", str(self.download_threads),
            "-t", str(self.extraction_threads),
            "--hide-download-progress",
        ]
        return cmd

    def download_one(self, srr: str) -> Optional[Path]:
        """Download a single SRR accession.

        Returns path to the output file on success, None on failure.
        """
        cmd = self._build_cmd([srr])
        logger.info(f"Downloading {srr}: {' '.join(cmd)}")

        ret = subprocess.run(cmd, capture_output=True, text=True)

        if ret.returncode != 0:
            logger.error(f"kingfisher failed for {srr}: {ret.stderr.strip()}")
            return None

        # Find output file(s)
        ext = self.output_format.replace(".", "*.")
        matches = sorted(self.output_dir.glob(f"{srr}*.{ext}"))
        if not matches:
            # Try broader match
            matches = sorted(self.output_dir.glob(f"{srr}*"))
            matches = [m for m in matches if m.is_file()]

        if matches:
            logger.info(f"  -> {[f.name for f in matches]}")
            return matches[0]

        logger.warning(f"No output files found for {srr}")
        return None

    def download_batch(
        self,
        metadata: pd.DataFrame,
        accession_col: str = "srr_accession",
        limit: int = 0,
        batch_size: int = 1,
    ) -> pd.DataFrame:
        """Download runs listed in a metadata DataFrame.

        Args:
            metadata: DataFrame with at least an accession column.
            accession_col: Column name containing SRR accessions.
            limit: Max number of runs to download (0 = all).
            batch_size: Number of accessions per kingfisher call.
                        Use 1 for per-run tracking, or higher for throughput.

        Returns:
            Copy of metadata with added 'fastq_path' and 'download_status' columns.
        """
        df = metadata.copy()
        if limit > 0:
            df = df.head(limit)

        paths = []
        statuses = []

        accessions = list(df[accession_col])

        for i in range(0, len(accessions), batch_size):
            batch = accessions[i : i + batch_size]
            logger.info(
                f"[{i + 1}–{i + len(batch)}/{len(accessions)}] "
                f"Downloading {batch}"
            )

            if batch_size == 1:
                result = self.download_one(batch[0])
                paths.append(str(result) if result else "")
                statuses.append("ok" if result else "failed")
            else:
                # Batch call — kingfisher handles multiple accessions
                cmd = self._build_cmd(batch)
                ret = subprocess.run(cmd, capture_output=True, text=True)
                for srr in batch:
                    matches = sorted(self.output_dir.glob(f"{srr}*"))
                    matches = [m for m in matches if m.is_file()]
                    if ret.returncode == 0 and matches:
                        paths.append(str(matches[0]))
                        statuses.append("ok")
                    else:
                        paths.append("")
                        statuses.append("failed")

        df["fastq_path"] = paths
        df["download_status"] = statuses

        ok = sum(1 for s in statuses if s == "ok")
        logger.info(f"Download complete: {ok}/{len(df)} succeeded")
        return df

    def download_from_list(
        self,
        accession_file: Path,
    ) -> subprocess.CompletedProcess:
        """Download all accessions listed in a text file (one per line).

        This uses kingfisher's native --run-identifiers-list for maximum
        efficiency on large batches.
        """
        cmd = [
            "kingfisher", "get",
            "--run-identifiers-list", str(accession_file),
            "-m", *self.methods,
            "-f", self.output_format,
            "--output-directory", str(self.output_dir),
            "--download-threads", str(self.download_threads),
            "-t", str(self.extraction_threads),
        ]
        logger.info(f"Batch download from {accession_file}: {' '.join(cmd)}")
        return subprocess.run(cmd, text=True)
