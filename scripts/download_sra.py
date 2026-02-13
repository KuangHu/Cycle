#!/usr/bin/env python3
"""Batch download SRA runs using kingfisher, from a metadata TSV or accession list."""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from cycle.download_manager import SRADownloader
from cycle.download_manager.config import DEFAULT_DOWNLOAD_METHODS


def main():
    parser = argparse.ArgumentParser(
        description="Download SRA runs using kingfisher."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--metadata",
        help="Path to metadata TSV (output of search_sra.py).",
    )
    input_group.add_argument(
        "--accession-list",
        help="Text file with one SRR accession per line.",
    )

    parser.add_argument(
        "-o", "--output-dir",
        default="data/sra_downloads",
        help="Download directory. Default: data/sra_downloads",
    )
    parser.add_argument(
        "-m", "--methods",
        nargs="+",
        default=DEFAULT_DOWNLOAD_METHODS,
        choices=["ena-ftp", "ena-ascp", "aws-http", "aws-cp", "gcp-cp", "prefetch"],
        help=f"Download methods in priority order. Default: {DEFAULT_DOWNLOAD_METHODS}",
    )
    parser.add_argument(
        "-f", "--format",
        default="fastq.gz",
        choices=["fastq", "fastq.gz", "fasta", "fasta.gz", "sra"],
        help="Output format. Default: fastq.gz",
    )
    parser.add_argument(
        "--download-threads",
        type=int, default=8,
        help="Connection threads for downloading. Default: 8",
    )
    parser.add_argument(
        "--extraction-threads",
        type=int, default=8,
        help="Threads for SRA extraction. Default: 8",
    )
    parser.add_argument(
        "--limit",
        type=int, default=0,
        help="Max runs to download (0 = all). Default: 0",
    )
    parser.add_argument(
        "--batch-size",
        type=int, default=1,
        help="Accessions per kingfisher call. Default: 1",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    downloader = SRADownloader(
        output_dir=args.output_dir,
        methods=args.methods,
        output_format=args.format,
        download_threads=args.download_threads,
        extraction_threads=args.extraction_threads,
    )

    if args.accession_list:
        # Direct kingfisher batch mode using accession list file
        downloader.download_from_list(Path(args.accession_list))
    else:
        # DataFrame-tracked mode from metadata TSV
        df = pd.read_csv(args.metadata, sep="\t")
        print(f"Loaded {len(df)} runs from {args.metadata}")

        result = downloader.download_batch(
            df, limit=args.limit, batch_size=args.batch_size,
        )

        out_tsv = str(Path(args.output_dir) / "download_status.tsv")
        result.to_csv(out_tsv, sep="\t", index=False)
        print(f"\nDownload status saved to {out_tsv}")


if __name__ == "__main__":
    main()
