#!/usr/bin/env python3
"""Detect partial IS circular intermediates via split-read back-jump analysis.

Maps reads to single-copy IS references and detects supplementary alignments
that jump backward on the reference — the signature of a sub-region of the IS
element circularizing.

Example
-------
python scripts/run_partial_circle.py \
    --is-records /path/to/sample_is_records_guide.json \
    --fastq /path/to/sample.fastq.gz \
    --sample-id SRR123456 \
    --output-dir /path/to/partial_circle_output \
    --threads 8
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.partial_circle import PartialCircleDetector


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--is-records",
        required=True,
        help="Path to *_is_records_guide.json with IS element sequences",
    )
    parser.add_argument(
        "--fastq",
        required=True,
        help="Path to sample FASTQ file (.fastq or .fastq.gz)",
    )
    parser.add_argument(
        "--sample-id",
        required=True,
        help="Sample accession (e.g. SRR123456)",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for partial circle output files",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=50,
        help="Min bp aligned on each side of junction (default: 50)",
    )
    parser.add_argument(
        "--min-circle-size",
        type=int,
        default=100,
        help="Min circle [S,E] span in bp (default: 100)",
    )
    parser.add_argument(
        "--max-circle-fraction",
        type=float,
        default=0.90,
        help="Max circle fraction of IS length (default: 0.90)",
    )
    parser.add_argument(
        "--breakpoint-tolerance",
        type=int,
        default=20,
        help="Breakpoint clustering tolerance in bp (default: 20)",
    )
    parser.add_argument(
        "--min-supporting-reads",
        type=int,
        default=2,
        help="Min supporting reads per call (default: 2)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Minimap2/samtools threads (default: 8)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    detector = PartialCircleDetector(
        output_dir=args.output_dir,
        min_overlap=args.min_overlap,
        min_circle_size=args.min_circle_size,
        max_circle_fraction=args.max_circle_fraction,
        breakpoint_tolerance=args.breakpoint_tolerance,
        min_supporting_reads=args.min_supporting_reads,
        threads=args.threads,
    )

    result = detector.run_sample(
        is_records_path=Path(args.is_records),
        sample_id=args.sample_id,
        fastq_path=Path(args.fastq),
    )

    log = logging.getLogger(__name__)
    if result:
        log.info("Done — results at %s", result)
    else:
        log.warning("No output produced for %s", args.sample_id)
        sys.exit(1)


if __name__ == "__main__":
    main()
