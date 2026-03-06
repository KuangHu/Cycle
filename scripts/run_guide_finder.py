#!/usr/bin/env python3
"""Standalone CLI for guide RNA alignment search on annotated IS records.

Usage:
    python scripts/run_guide_finder.py \
        --input-dir /path/to/is_formatter_output/ \
        --parallel 20
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.guide_finder import GuideFinder

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Search for guide-like alignments between flanking and noncoding regions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing */*_is_records_annotated.json files",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel workers. Default: 1",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=9,
        help="Minimum alignment length. Default: 9",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=1,
        help="Maximum mismatches allowed. Default: 1",
    )
    parser.add_argument(
        "--sample-list",
        default=None,
        help="TSV file with 'srr_accession' column to restrict processing to those samples.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        logger.error(f"Input directory not found: {input_dir}")
        sys.exit(1)

    sample_ids = None
    if args.sample_list:
        import pandas as pd
        sample_ids = set(pd.read_csv(args.sample_list, sep="\t")["srr_accession"])
        logger.info(f"Restricting to {len(sample_ids)} samples from {args.sample_list}")

    gf = GuideFinder(
        min_length=args.min_length,
        max_mismatches=args.max_mismatches,
    )
    results = gf.find_guides_batch(input_dir, parallel=args.parallel, sample_ids=sample_ids)

    succeeded = sum(1 for v in results.values() if v is not None)
    logger.info(f"Done: {succeeded}/{len(results)} samples processed successfully")


if __name__ == "__main__":
    main()
