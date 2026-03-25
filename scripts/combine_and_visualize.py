#!/usr/bin/env python3
"""Combine batches 006-010, find systems with guide + circle evidence,
generate PNG + GBK for each.

Usage:
    python scripts/combine_and_visualize.py \
        --batches 6 7 8 9 10 \
        --output-dir /groups/rubin/projects/kuang/out/IS_cycle/combined_006_010 \
        --parallel 20
"""

import argparse
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cycle.visualizer.visualizer import ISElementVisualizer
from cycle.visualizer.genbank import ISElementGenBank


def load_partial_circles(batch_dirs):
    """Load all partial circle calls across batches, keyed by is_id."""
    pc_by_is = {}
    for bdir in batch_dirs:
        pc_dir = os.path.join(bdir, "partial_circle_output")
        for sj in glob.glob(
            os.path.join(pc_dir, "**/*_partial_circle_summary.json"),
            recursive=True,
        ):
            with open(sj) as f:
                data = json.load(f)
            if isinstance(data, list):
                for call in data:
                    is_id = call.get("is_id", "")
                    if is_id:
                        pc_by_is.setdefault(is_id, []).append(call)
    return pc_by_is


def collect_records(batch_dirs, min_orf_length=100):
    """Collect all IS records with guide + circle evidence.

    Filters out records where the longest ORF is shorter than *min_orf_length* bp.
    """
    records = []
    skipped_short_orf = 0
    for bdir in batch_dirs:
        batch_name = os.path.basename(bdir)
        fmt_dir = os.path.join(bdir, "is_formatter_output")
        for gj in sorted(
            glob.glob(os.path.join(fmt_dir, "*/*_is_records_guide.json"))
        ):
            with open(gj) as f:
                data = json.load(f)
            for rec in data:
                ce = rec.get("circle_evidence", {})
                gs = rec.get("guide_summary", {})
                n_th = ce.get("n_tail_head_reads", 0)
                n_guide = gs.get("n_hits", 0)
                if n_guide > 0 and n_th > 0:
                    # Filter by longest ORF length
                    orfs = (rec.get("orf_annotation") or {}).get("orfs", [])
                    max_orf_len = max((o.get("length_nt", 0) for o in orfs), default=0)
                    if max_orf_len < min_orf_length:
                        skipped_short_orf += 1
                        continue
                    rec["_batch"] = batch_name
                    records.append(rec)
    if skipped_short_orf:
        print(f"  Skipped {skipped_short_orf} records with longest ORF < {min_orf_length}bp")
    return records


def visualize_one(rec, pc_by_is, output_dir, dpi=150):
    """Generate PNG + GBK for one record."""
    is_id = rec["is_id"]
    sample_id = rec["sample_id"]
    tag = f"{sample_id}__{is_id[:8]}"

    png_path = os.path.join(output_dir, f"{tag}.png")
    gbk_path = os.path.join(output_dir, f"{tag}.gbk")

    if os.path.exists(png_path) and os.path.exists(gbk_path):
        return tag, True

    alignments = rec.get("guide_hits", [])
    circle_info = rec.get("circle_evidence")
    partial_circles = pc_by_is.get(is_id, None)

    try:
        vis = ISElementVisualizer()
        vis.save_element_png(
            rec, alignments, png_path,
            dpi=dpi, circle_info=circle_info,
            partial_circles=partial_circles,
        )

        gb = ISElementGenBank()
        gb.save_genbank(
            rec, alignments, gbk_path,
            circle_info=circle_info,
            partial_circles=partial_circles,
        )
        return tag, True
    except Exception as e:
        return tag, f"ERROR: {e}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", nargs="+", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--parallel", type=int, default=20)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--min-orf-length", type=int, default=100,
                        help="Minimum longest ORF length in bp (default: 100)")
    args = parser.parse_args()

    base = "/groups/rubin/projects/kuang/out/IS_cycle"
    batch_dirs = [os.path.join(base, f"batch_{b:03d}") for b in args.batches]

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading partial circle data...")
    pc_by_is = load_partial_circles(batch_dirs)
    print(f"  {len(pc_by_is)} IS elements with partial circle calls")

    print("Collecting records with guide + circle evidence...")
    records = collect_records(batch_dirs, min_orf_length=args.min_orf_length)
    print(f"  {len(records)} systems to visualize")

    # Save summary TSV
    summary_path = os.path.join(args.output_dir, "systems_summary.tsv")
    with open(summary_path, "w") as f:
        f.write("sample_id\tis_id\tbatch\tn_th_reads\tn_gh_reads\tn_tg_reads\t"
                "n_guide_hits\thas_partial_circle\n")
        for rec in records:
            ce = rec.get("circle_evidence", {})
            gs = rec.get("guide_summary", {})
            is_id = rec["is_id"]
            f.write(f"{rec['sample_id']}\t{is_id}\t{rec.get('_batch','')}\t"
                    f"{ce.get('n_tail_head_reads',0)}\t"
                    f"{ce.get('n_genome_head_reads',0)}\t"
                    f"{ce.get('n_tail_genome_reads',0)}\t"
                    f"{gs.get('n_hits',0)}\t"
                    f"{is_id in pc_by_is}\n")
    print(f"  Summary: {summary_path}")

    print(f"Generating PNGs + GBKs ({args.parallel} parallel)...")
    done = 0
    errors = 0

    if args.parallel <= 1:
        for rec in records:
            tag, result = visualize_one(rec, pc_by_is, args.output_dir, args.dpi)
            done += 1
            if result is not True:
                errors += 1
                print(f"  {tag}: {result}")
            if done % 100 == 0:
                print(f"  {done}/{len(records)} done")
    else:
        with ProcessPoolExecutor(max_workers=args.parallel) as pool:
            futures = {
                pool.submit(visualize_one, rec, pc_by_is, args.output_dir, args.dpi): rec
                for rec in records
            }
            for fut in as_completed(futures):
                tag, result = fut.result()
                done += 1
                if result is not True:
                    errors += 1
                    print(f"  {tag}: {result}")
                if done % 500 == 0:
                    print(f"  {done}/{len(records)} done")

    print(f"\nDone: {done} systems, {errors} errors")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
