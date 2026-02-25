#!/usr/bin/env python
"""
Find IS elements with IS110 transposase AND tail-head circular junctions.

Usage:
    python scripts/run_is110_filter.py \
        --protein-fasta /path/to/all_proteins.faa \
        --input-dir /path/to/is_formatter_output \
        --output-dir /path/to/is110_circular \
        --cpus 8 \
        --min-tail-head 1
"""

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cycle.is110_filter import IS110Filter


def main():
    parser = argparse.ArgumentParser(
        description="Filter IS elements for IS110 protein + tail-head circular junctions.",
    )
    parser.add_argument(
        "--protein-fasta", required=True,
        help="Path to all_proteins.faa from system clustering.",
    )
    parser.add_argument(
        "--input-dir", required=True,
        help="Directory containing sample subdirs with *_is_records_guide.json.",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for filtered results.",
    )
    parser.add_argument(
        "--cpus", type=int, default=8,
        help="Threads for hmmsearch (default: 8).",
    )
    parser.add_argument(
        "--min-tail-head", type=int, default=1,
        help="Minimum n_tail_head_reads to pass filter (default: 1).",
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="hmmsearch E-value threshold (default: 1e-5).",
    )
    parser.add_argument(
        "--skip-hmmsearch", action="store_true",
        help="Skip hmmsearch; reuse existing DEDD_hits.tbl and Tnp20_hits.tbl in output-dir.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    filt = IS110Filter(evalue=args.evalue)

    if args.skip_hmmsearch:
        logging.info("Reusing existing hmmsearch results in %s", args.output_dir)
        is110_ids = set()
        for name in ("DEDD", "Tnp20"):
            tbl = os.path.join(args.output_dir, f"{name}_hits.tbl")
            if os.path.exists(tbl):
                ids = filt.parse_tblout(tbl)
                logging.info("%s: %d hits → %d transposons", name, filt._last_n_hits, len(ids))
                is110_ids.update(ids)
        logging.info("Combined: %d IS110 transposons", len(is110_ids))
        records = filt.filter_records(args.input_dir, is110_ids, min_tail_head=args.min_tail_head)
        filt.export_results(records, args.output_dir)
    else:
        filt.run(
            protein_fasta=args.protein_fasta,
            formatter_dir=args.input_dir,
            output_dir=args.output_dir,
            cpus=args.cpus,
            min_tail_head=args.min_tail_head,
        )

    logging.info("Done.")


if __name__ == "__main__":
    main()
