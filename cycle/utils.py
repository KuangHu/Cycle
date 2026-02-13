"""Shared utility functions used across Cycle subpackages."""

import re
from pathlib import Path


def slugify(name: str) -> str:
    """Convert organism name to a filesystem-safe slug."""
    slug = name.lower().strip()
    slug = re.sub(r"[^a-z0-9]+", "_", slug)
    return slug.strip("_")


def find_fastq(fastq_dir: Path, sample_id: str) -> Path | None:
    """Find FASTQ file for a sample, trying common naming patterns."""
    for pattern in (
        f"{sample_id}.fastq.gz",
        f"{sample_id}.fastq",
        f"{sample_id}_1.fastq.gz",
        f"{sample_id}_pass.fastq.gz",
    ):
        path = fastq_dir / pattern
        if path.exists():
            return path

    matches = sorted(fastq_dir.glob(f"{sample_id}*fastq*"))
    if matches:
        return matches[0]

    return None
