#!/usr/bin/env python3
"""Cluster transposon systems by shared protein families.

Reads *_is_records_guide.json files from formatter output directories,
clusters proteins with MMseqs2, groups transposons by shared protein
families (Louvain community detection), and classifies within-cluster
variants by flanking region and guide sequence.

Example
-------
python scripts/run_system_clustering.py \
    --input-dirs /path/to/batch_000/is_formatter_output \
                 /path/to/batch_001/is_formatter_output \
    --output-dir /path/to/system_clustering \
    --threads 48
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.system_clustering import SystemClusterer


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
        "--output-dir",
        required=True,
        help="Directory for clustering output files",
    )
    parser.add_argument(
        "--min-seq-id",
        type=float,
        default=0.3,
        help="MMseqs2 minimum sequence identity (default: 0.3)",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=0.8,
        help="MMseqs2 coverage threshold (default: 0.8)",
    )
    parser.add_argument(
        "--louvain-resolution",
        type=float,
        default=1.0,
        help="Louvain resolution parameter (default: 1.0)",
    )
    parser.add_argument(
        "--jaccard-sim-threshold",
        type=float,
        default=0.3,
        help="Minimum Jaccard similarity for graph edges (default: 0.3)",
    )
    parser.add_argument(
        "--flanking-edit-threshold",
        type=int,
        default=10,
        help="Max edit distance for same L1 flanking group (default: 10)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="MMseqs2 threads (default: 8)",
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

    clusterer = SystemClusterer(
        output_dir=args.output_dir,
        min_seq_id=args.min_seq_id,
        coverage=args.coverage,
        louvain_resolution=args.louvain_resolution,
        jaccard_sim_threshold=args.jaccard_sim_threshold,
        flanking_edit_threshold=args.flanking_edit_threshold,
        mmseqs_threads=args.threads,
    )

    results_path = clusterer.run(
        formatter_dirs=args.input_dirs,
        skip_existing=not args.no_skip,
    )
    logging.getLogger(__name__).info("Done — results at %s", results_path)


if __name__ == "__main__":
    main()
