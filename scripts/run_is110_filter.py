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
        "--max-orf-gap", type=int, default=300,
        help="Max bp gap between adjacent ORFs for two-protein IS110 (default: 300).",
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

    filt = IS110Filter(evalue=args.evalue, max_orf_gap=args.max_orf_gap)

    if args.skip_hmmsearch:
        logging.info("Reusing existing hmmsearch results in %s", args.output_dir)
        # Parse both tblout files, then apply the same two-case logic
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

        # Apply same two-case logic as run_hmmsearch
        from collections import defaultdict
        trans_domains = defaultdict(lambda: {"DEDD": [], "Tnp20": []})
        for domain, hits in protein_hits.items():
            for prot_id, (trans_id, start, end, strand) in hits.items():
                trans_domains[trans_id][domain].append((prot_id, start, end, strand))

        is110_ids = set()
        n_case1 = n_case2 = 0
        for trans_id, domains in trans_domains.items():
            dedd_prots = {p[0] for p in domains["DEDD"]}
            tnp20_prots = {p[0] for p in domains["Tnp20"]}
            if dedd_prots & tnp20_prots:
                is110_ids.add(trans_id)
                n_case1 += 1
                continue
            if dedd_prots and tnp20_prots:
                if filt._has_adjacent_pair(domains["DEDD"], domains["Tnp20"]):
                    is110_ids.add(trans_id)
                    n_case2 += 1

        logging.info(
            "IS110: %d transposons (case1=%d same-ORF, case2=%d adjacent-ORFs, "
            "rejected %d single-domain)",
            len(is110_ids), n_case1, n_case2,
            len(trans_domains) - len(is110_ids),
        )
        records = filt.filter_records(
            args.input_dir, is110_ids,
            min_tail_head=args.min_tail_head,
            trans_domains=dict(trans_domains),
        )
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
