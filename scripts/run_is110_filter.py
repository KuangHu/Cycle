#!/usr/bin/env python
"""
Find IS elements with IS110 transposase, split by circle evidence.

Usage:
    python scripts/run_is110_filter.py \
        --protein-fasta /path/to/all_proteins.faa \
        --input-dirs /path/to/batch_006/is_formatter_output \
                     /path/to/batch_007/is_formatter_output \
        --output-dir /path/to/is110_output \
        --partial-circle-dirs /path/to/batch_006/partial_circle_output \
                              /path/to/batch_007/partial_circle_output \
        --cpus 8
"""

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cycle.is110_filter import IS110Filter


def main():
    parser = argparse.ArgumentParser(
        description="Filter IS elements for IS110 protein, split by circle evidence.",
    )
    parser.add_argument(
        "--protein-fasta", required=True,
        help="Path to all_proteins.faa.",
    )
    parser.add_argument(
        "--input-dirs", nargs="+", required=True,
        help="Formatter output directories containing */*_is_records_guide.json.",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for results.",
    )
    parser.add_argument(
        "--partial-circle-dirs", nargs="+", default=None,
        help="Directories with *_partial_circle_summary.json files.",
    )
    parser.add_argument(
        "--cpus", type=int, default=8,
        help="Threads for hmmsearch (default: 8).",
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="hmmsearch E-value threshold (default: 1e-5).",
    )
    parser.add_argument(
        "--max-orf-gap", type=int, default=300,
        help="Max bp gap between adjacent ORFs for two-protein IS110 (default: 300).",
    )
    parser.add_argument(
        "--skip-hmmsearch", action="store_true",
        help="Skip hmmsearch; reuse existing DEDD_hits.tbl and Tnp20_hits.tbl in output-dir.",
    )
    parser.add_argument(
        "--no-visualize", action="store_true",
        help="Skip PNG + GBK generation.",
    )
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="PNG resolution (default: 150).",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    filt = IS110Filter(evalue=args.evalue, max_orf_gap=args.max_orf_gap)

    if args.skip_hmmsearch:
        logging.info("Reusing existing hmmsearch results in %s", args.output_dir)
        protein_hits = {}
        for name in ("DEDD", "Tnp20"):
            tbl = os.path.join(args.output_dir, f"{name}_hits.tbl")
            if os.path.exists(tbl):
                hits = filt.parse_tblout(tbl)
                protein_hits[name] = hits
                logging.info("%s: %d protein hits", name, len(hits))
            else:
                protein_hits[name] = {}
                logging.warning("%s tblout not found: %s", name, tbl)

        is110_ids, trans_domains = filt._identify_is110(protein_hits)

        # Load partial circle data
        pc_by_is = {}
        partial_circle_ids = set()
        if args.partial_circle_dirs:
            pc_by_is = filt.load_partial_circles(args.partial_circle_dirs)
            partial_circle_ids = set(pc_by_is.keys())

        # Filter and split
        all_with = []
        all_without = []
        for fmt_dir in args.input_dirs:
            w, wo = filt.filter_records(
                fmt_dir, is110_ids,
                trans_domains=trans_domains,
                partial_circle_ids=partial_circle_ids,
            )
            all_with.extend(w)
            all_without.extend(wo)

        logging.info("Total: %d with circle, %d without", len(all_with), len(all_without))
        filt.export_results(all_with, all_without, args.output_dir)

        if not args.no_visualize:
            wc_dir = os.path.join(args.output_dir, "with_circle_evidence")
            wo_dir = os.path.join(args.output_dir, "without_circle_evidence")
            filt.visualize(all_with, wc_dir, pc_by_is=pc_by_is, dpi=args.dpi)
            filt.visualize(all_without, wo_dir, pc_by_is=pc_by_is, dpi=args.dpi)
    else:
        filt.run(
            protein_fasta=args.protein_fasta,
            formatter_dirs=args.input_dirs,
            output_dir=args.output_dir,
            cpus=args.cpus,
            partial_circle_dirs=args.partial_circle_dirs,
            visualize=not args.no_visualize,
            dpi=args.dpi,
        )

    logging.info("Done.")


if __name__ == "__main__":
    main()
