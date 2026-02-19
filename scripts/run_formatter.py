#!/usr/bin/env python3
"""Standalone test script for IS element formatting.

Usage:
    python scripts/run_formatter.py \
      --circle-dir /path/to/circle_output/organism_slug \
      --organism "Organism name" \
      --outdir /tmp/test_format \
      --flank-size 80
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.is_formatter import ISFormatter
from cycle.is_formatter.config import (
    DEFAULT_ASSEMBLY_TIMEOUT,
    DEFAULT_FLANK_SIZE,
    DEFAULT_FORMATTER_OUTPUT_DIR,
    DEFAULT_MIN_READS_FOR_ASSEMBLY,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract IS elements + flanking regions from circle detection output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--circle-dir", required=True,
        help="Path to the organism's circle detection output directory "
             "(contains *_circle_summary.tsv and *.circle.sorted.bam).",
    )
    parser.add_argument(
        "--organism", required=True,
        help="Organism name (used for output naming).",
    )
    parser.add_argument(
        "--outdir", default=DEFAULT_FORMATTER_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_FORMATTER_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--flank-size", type=int, default=DEFAULT_FLANK_SIZE,
        help=f"Flanking region size in bp. Default: {DEFAULT_FLANK_SIZE}",
    )
    parser.add_argument(
        "--min-reads", type=int, default=DEFAULT_MIN_READS_FOR_ASSEMBLY,
        help=f"Minimum reads for assembly. Default: {DEFAULT_MIN_READS_FOR_ASSEMBLY}",
    )
    parser.add_argument(
        "--assembly-timeout", type=int, default=DEFAULT_ASSEMBLY_TIMEOUT,
        help=f"Timeout (seconds) per IS element assembly. Default: {DEFAULT_ASSEMBLY_TIMEOUT}",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Threads for minimap2/minipolish. Default: 4",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    circle_dir = Path(args.circle_dir)
    if not circle_dir.exists():
        logging.error(f"Circle directory does not exist: {circle_dir}")
        sys.exit(1)

    formatter = ISFormatter(
        output_dir=args.outdir,
        flank_size=args.flank_size,
        min_reads=args.min_reads,
        assembly_timeout=args.assembly_timeout,
        threads=args.threads,
    )

    result = formatter.run_organism(
        circle_dir=circle_dir,
        organism=args.organism,
    )

    if result:
        logging.info(f"Output written to: {result}")
    else:
        logging.warning("No output produced.")
        sys.exit(1)


if __name__ == "__main__":
    main()
