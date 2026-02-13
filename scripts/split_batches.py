#!/usr/bin/env python3
"""Split master SRA metadata into diverse batches of ~1000 samples.

Strategy: round-robin assignment by organism so that each batch contains
samples from as many different species as possible, with both platforms
(Nanopore / PacBio) represented proportionally.

Usage:
    python scripts/split_batches.py \
        --metadata data/sra_metadata.tsv \
        --output-dir data/batches \
        --batch-size 1000
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def split_batches(
    metadata: pd.DataFrame,
    batch_size: int = 1000,
) -> list[pd.DataFrame]:
    """Split metadata into diverse batches via round-robin by organism.

    Within each organism, samples are shuffled so that platform and other
    attributes are distributed randomly.  Organisms are then interleaved
    round-robin across batches, producing batches where each one samples
    broadly across the full species range.
    """
    # Shuffle within each organism (mixes platforms, projects, etc.)
    metadata = metadata.sample(frac=1, random_state=42).reset_index(drop=True)

    # Sort organisms by count (descending) so large groups spread evenly
    org_counts = metadata["organism"].value_counts()
    org_order = org_counts.index.tolist()

    n_batches = (len(metadata) + batch_size - 1) // batch_size

    # Assign batch index round-robin within each organism
    batch_assignments = pd.Series(index=metadata.index, dtype=int)
    batch_cursors = [0] * n_batches  # track how full each batch is

    for org in org_order:
        org_mask = metadata["organism"] == org
        org_indices = metadata.index[org_mask].tolist()

        # Find the least-full batches and assign round-robin
        for i, idx in enumerate(org_indices):
            # Pick the batch with fewest samples so far
            target = min(range(n_batches), key=lambda b: batch_cursors[b])
            batch_assignments[idx] = target
            batch_cursors[target] += 1

    # Build batch DataFrames
    batches = []
    for b in range(n_batches):
        batch_df = metadata.loc[batch_assignments == b].copy()
        batch_df = batch_df.sample(frac=1, random_state=b).reset_index(drop=True)
        batches.append(batch_df)

    return batches


def main():
    parser = argparse.ArgumentParser(
        description="Split SRA metadata into diverse batches.",
    )
    parser.add_argument(
        "--metadata", default="data/sra_metadata.tsv",
        help="Path to master metadata TSV.",
    )
    parser.add_argument(
        "--output-dir", default="data/batches",
        help="Output directory for batch TSVs.",
    )
    parser.add_argument(
        "--batch-size", type=int, default=1000,
        help="Target samples per batch. Default: 1000",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    metadata = pd.read_csv(args.metadata, sep="\t")
    logger.info(
        f"Loaded {len(metadata)} samples, "
        f"{metadata['organism'].nunique()} organisms"
    )

    batches = split_batches(metadata, batch_size=args.batch_size)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for i, batch_df in enumerate(batches):
        batch_path = out_dir / f"batch_{i:03d}.tsv"
        batch_df.to_csv(batch_path, sep="\t", index=False)

    # Print summary
    logger.info(f"Wrote {len(batches)} batches to {out_dir}/")
    for i, batch_df in enumerate(batches):
        n_org = batch_df["organism"].nunique()
        platforms = batch_df["platform"].value_counts().to_dict()
        plat_str = ", ".join(f"{k}={v}" for k, v in platforms.items())
        logger.info(f"  batch_{i:03d}: {len(batch_df)} samples, {n_org} organisms ({plat_str})")


if __name__ == "__main__":
    main()
