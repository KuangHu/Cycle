#!/usr/bin/env python
"""
CLI script to generate PNG diagrams and/or GenBank files for IS elements
from *_is_records_guide.json files produced by the Cycle pipeline.

Usage:
    python scripts/run_visualizer.py \
        --input-dir /path/to/is_formatter_output \
        --output-dir /path/to/visualizations \
        --format both \
        --parallel 20 \
        --dpi 150
"""

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cycle.visualizer import ISElementVisualizer, ISElementGenBank


def main():
    parser = argparse.ArgumentParser(
        description="Generate PNG diagrams and GenBank files for IS elements.",
    )
    parser.add_argument(
        "--input-dir", required=True,
        help="Directory containing sample subdirs with *_is_records_guide.json files.",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory (default: <input-dir>/../visualizations).",
    )
    parser.add_argument(
        "--format", choices=["png", "gbk", "both"], default="both",
        help="Output format (default: both).",
    )
    parser.add_argument(
        "--parallel", type=int, default=1,
        help="Number of parallel workers (default: 1).",
    )
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="PNG resolution (default: 150).",
    )
    parser.add_argument(
        "--figure-width", type=int, default=12,
        help="PNG figure width in inches (default: 12).",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    input_dir = os.path.abspath(args.input_dir)
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
    else:
        output_dir = os.path.join(os.path.dirname(input_dir), "visualizations")

    fmt = args.format

    if fmt in ("png", "both"):
        logging.info("Generating PNG diagrams → %s", output_dir)
        vis = ISElementVisualizer()
        vis.visualize_batch(
            input_dir, output_dir,
            parallel=args.parallel,
            dpi=args.dpi,
            figure_width=args.figure_width,
        )

    if fmt in ("gbk", "both"):
        logging.info("Generating GenBank files → %s", output_dir)
        gb = ISElementGenBank()
        gb.export_batch(
            input_dir, output_dir,
            parallel=args.parallel,
        )

    logging.info("Done.")


if __name__ == "__main__":
    main()
