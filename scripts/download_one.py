#!/usr/bin/env python3
"""Download a single SRA accession by array index.

Designed for SLURM job arrays:
    sbatch --array=0-499 ... --wrap "python scripts/download_one.py --metadata batch.tsv --index \$SLURM_ARRAY_TASK_ID --outdir /path/to/downloads"
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from cycle.download_manager import SRADownloader


def main():
    parser = argparse.ArgumentParser(
        description="Download a single SRA accession by row index.",
    )
    parser.add_argument("--metadata", required=True, help="Metadata TSV file")
    parser.add_argument("--index", type=int, required=True, help="Row index (from SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--outdir", required=True, help="Output directory for FASTQ files")
    parser.add_argument("--accession-col", default="srr_accession", help="Column with accession IDs")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    df = pd.read_csv(args.metadata, sep="\t")
    if args.index >= len(df):
        logging.info(f"Index {args.index} >= {len(df)} rows, nothing to do")
        sys.exit(0)

    acc = df.iloc[args.index][args.accession_col]

    # Skip if already downloaded
    outdir = Path(args.outdir)
    if list(outdir.glob(f"{acc}*fastq*")):
        logging.info(f"{acc}: already exists, skipping")
        sys.exit(0)

    dl = SRADownloader(output_dir=args.outdir, download_threads=1, extraction_threads=1)
    result = dl.download_one(acc)

    if result:
        logging.info(f"{acc}: OK -> {result}")
    else:
        logging.error(f"{acc}: FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
