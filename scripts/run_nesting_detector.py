#!/usr/bin/env python3
"""Detect IS-within-IS nesting events among IS110 elements.

Runs pairwise alignment of IS110 extended sequences (flanking + IS) to find
elements that are nearly identical except one has extra sequence inserted —
evidence of IS-within-IS nesting.

Example
-------
    python scripts/run_nesting_detector.py \
        --records is110_circular_batch_000/is110_circular_records.json \
        --output-dir nesting_output
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.nesting_detector import NestingDetector


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--records",
        required=True,
        help="Path to is110_circular_records.json",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for nesting detection output",
    )
    parser.add_argument(
        "--mm2-preset",
        default="asm10",
        help="minimap2 preset (default: asm10)",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.90,
        help="Minimum alignment block identity (default: 0.90)",
    )
    parser.add_argument(
        "--min-insertion-size",
        type=int,
        default=50,
        help="Minimum insertion size in bp (default: 50)",
    )
    parser.add_argument(
        "--min-block-length",
        type=int,
        default=100,
        help="Minimum aligned block length in bp (default: 100)",
    )
    parser.add_argument(
        "--flanking-pad",
        type=int,
        default=80,
        help="Flanking bp included in extended sequence (default: 80)",
    )
    parser.add_argument(
        "--min-length-ratio",
        type=float,
        default=1.02,
        help="Host must be >= this ratio longer than core (default: 1.02)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="minimap2 threads (default: 4)",
    )
    parser.add_argument(
        "--no-skip",
        action="store_true",
        help="Re-run even if output already exists",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    detector = NestingDetector(
        output_dir=args.output_dir,
        mm2_preset=args.mm2_preset,
        min_identity=args.min_identity,
        min_insertion_size=args.min_insertion_size,
        min_block_length=args.min_block_length,
        flanking_pad=args.flanking_pad,
        min_length_ratio=args.min_length_ratio,
        threads=args.threads,
    )

    detector.run(
        records_path=Path(args.records),
        skip_existing=not args.no_skip,
    )

    logging.getLogger(__name__).info("Done — results at %s", args.output_dir)


if __name__ == "__main__":
    main()
