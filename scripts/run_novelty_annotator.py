#!/usr/bin/env python3
"""Annotate IS elements with ISfinder BLAST hits and score cluster novelty.

BLASTs all detected IS elements against the ISfinder reference database,
annotates each record with the closest known match, then computes a
composite novelty score per system cluster based on divergence from
ISfinder, within-cluster diversity, and mosaic structure.

Example
-------
python scripts/run_novelty_annotator.py \
    --input-dirs /path/to/batch_000/is_formatter_output \
                 /path/to/batch_001/is_formatter_output \
    --clusters /path/to/system_clusters.json \
    --output-dir /path/to/novelty_output \
    --threads 8
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.novelty_annotator import NoveltyAnnotator


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dirs",
        nargs="+",
        required=True,
        help="Formatter output directories containing */*_is_records_guide.json",
    )
    parser.add_argument(
        "--clusters",
        required=True,
        help="Path to system_clusters.json from system clustering",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for novelty annotation output files",
    )
    parser.add_argument(
        "--isfinder-fasta",
        default=None,
        help=(
            "Path to ISfinder FASTA (ISfinder_raw.fna). "
            "Auto-detected from is_reference/ sibling dir if omitted."
        ),
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="BLAST E-value threshold (default: 1e-5)",
    )
    parser.add_argument(
        "--max-target-seqs",
        type=int,
        default=5,
        help="Max BLAST target sequences per query (default: 5)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="BLAST threads (default: 8)",
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

    annotator = NoveltyAnnotator(
        output_dir=args.output_dir,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        threads=args.threads,
    )

    annotator.run(
        formatter_dirs=args.input_dirs,
        clusters_path=args.clusters,
        isfinder_fasta=args.isfinder_fasta,
        skip_existing=not args.no_skip,
    )

    logging.getLogger(__name__).info("Done — results at %s", args.output_dir)


if __name__ == "__main__":
    main()
