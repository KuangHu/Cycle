#!/usr/bin/env python3
"""Search NCBI SRA for bacterial nanopore WGS runs and save metadata."""

import argparse
import logging
import sys

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))

from cycle.download_manager import SRASearcher
from cycle.download_manager.config import DEFAULT_SRA_QUERY


def main():
    parser = argparse.ArgumentParser(
        description="Search NCBI SRA and collect run metadata."
    )
    parser.add_argument(
        "--email", required=True, help="Email for NCBI Entrez (required by NCBI)."
    )
    parser.add_argument("--api-key", default=None, help="NCBI API key (optional).")
    parser.add_argument(
        "--query",
        default=DEFAULT_SRA_QUERY,
        help=f"Entrez search query. Default: {DEFAULT_SRA_QUERY}",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=0,
        help="Max number of UIDs to retrieve (0 = all). Default: 0",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="data/sra_metadata.tsv",
        help="Output TSV path. Default: data/sra_metadata.tsv",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    searcher = SRASearcher(email=args.email, api_key=args.api_key)
    df = searcher.search_and_collect(
        query=args.query,
        max_results=args.max_results,
        output=args.output,
    )
    print(f"\nCollected {len(df)} runs -> {args.output}")
    print(df[["srr_accession", "organism", "total_bases"]].head(10).to_string())


if __name__ == "__main__":
    main()
