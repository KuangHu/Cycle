#!/usr/bin/env python3
"""Standalone CLI for ORF annotation of IS formatter output.

Usage:
    python scripts/run_orf_annotator.py \
        --input-dir /path/to/is_formatter_th_output/ \
        --parallel 4
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.orf_annotator import ORFAnnotator

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate IS element records with ORF predictions (Prodigal).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Formatter output directory containing *_is_records.json files",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel workers. Default: 1",
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

    annotator = ORFAnnotator()
    results = annotator.annotate_batch(input_dir, parallel=args.parallel)

    succeeded = sum(1 for v in results.values() if v is not None)
    logger.info(f"Done: {succeeded}/{len(results)} samples annotated successfully")


if __name__ == "__main__":
    main()
