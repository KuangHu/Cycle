#!/usr/bin/env python3
"""Check pipeline status across all batches.

Usage:
    python scripts/check_status.py                    # all batches
    python scripts/check_status.py 0 4                # batches 000-004
    python scripts/check_status.py --outroot /path    # custom output root
"""

import argparse
import sys
from pathlib import Path

def count_lines(path, skip_header=True):
    """Count non-empty lines in a file."""
    if not path.exists():
        return 0
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    return max(0, len(lines) - (1 if skip_header else 0))

def count_by_status(tsv_path, status_col, ok_val):
    """Count rows where status_col == ok_val."""
    if not tsv_path.exists():
        return 0, 0
    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
        if status_col not in header:
            return 0, 0
        idx = header.index(status_col)
        total = 0
        ok = 0
        for line in f:
            if not line.strip():
                continue
            total += 1
            fields = line.strip().split("\t")
            if len(fields) > idx and fields[idx] == ok_val:
                ok += 1
    return ok, total

def count_sample_dirs(step_dir, suffix):
    """Count sample-keyed output files: {step_dir}/{SRR|ERR|DRR}*/{*}{suffix}"""
    if not step_dir.exists():
        return 0
    count = 0
    for d in step_dir.iterdir():
        if not d.is_dir():
            continue
        name = d.name
        if not (name.startswith("SRR") or name.startswith("ERR") or name.startswith("DRR")):
            continue
        if list(d.glob(f"*{suffix}")):
            count += 1
    return count

def count_organism_dirs(step_dir, suffix):
    """Count organism-keyed output files (dirs that are NOT SRR/ERR/DRR)."""
    if not step_dir.exists():
        return 0
    count = 0
    for d in step_dir.iterdir():
        if not d.is_dir():
            continue
        name = d.name
        if name.startswith("SRR") or name.startswith("ERR") or name.startswith("DRR"):
            continue
        if list(d.glob(f"*{suffix}")):
            count += 1
    return count

def check_batch(batch_dir, batch_name, metadata_dir):
    """Check all pipeline stages for one batch."""
    bd = Path(batch_dir)
    meta_path = Path(metadata_dir) / f"{batch_name}.tsv"
    meta_sniffles = bd / "metadata_for_sniffles.tsv"

    # Expected samples
    expected = count_lines(meta_path) if meta_path.exists() else 0
    meta_sniffles_n = count_lines(meta_sniffles) if meta_sniffles.exists() else 0

    # Stage 1: Download
    dl_ok, dl_total = count_by_status(
        bd / "sra_downloads" / "download_status.tsv", "download_status", "ok"
    )

    # Stage 2: Resolve
    resolve_status = bd / "reference_genomes" / "resolve_status.tsv"
    resolve_n = count_lines(resolve_status) if resolve_status.exists() else 0
    # Also count .fna files on disk
    ref_dir = bd / "reference_genomes"
    fna_count = len(list(ref_dir.glob("*/*_genomic.fna"))) if ref_dir.exists() else 0

    # Stage 4: Align
    align_ok, align_total = count_by_status(
        bd / "alignments" / "alignment_status.tsv", "status", "ok"
    )

    # Stage 6: tldr (sample-keyed)
    tldr_dir = bd / "tldr_output"
    tldr_sample = count_sample_dirs(tldr_dir, ".table.txt")
    tldr_org = count_organism_dirs(tldr_dir, ".table.txt")

    # Stage 6alt: Sniffles (sample-keyed)
    sniffles_dir = bd / "sniffles_output"
    sniffles_sample = count_sample_dirs(sniffles_dir, ".table.txt")
    sniffles_org = count_organism_dirs(sniffles_dir, ".table.txt")

    # Stage 7: Circle (sample-keyed)
    circle_dir = bd / "circle_output"
    circle_sample = count_sample_dirs(circle_dir, "_circle_summary.tsv")
    circle_org = count_organism_dirs(circle_dir, "_circle_summary.tsv")

    # Stage 8: Format (sample-keyed)
    fmt_dir = bd / "is_formatter_output"
    fmt_sample = count_sample_dirs(fmt_dir, "_is_records.json")
    fmt_org = count_organism_dirs(fmt_dir, "_is_records.json")

    # Also check tldr-based formatter if exists
    fmt_tldr_dir = bd / "is_formatter_output_tldr"
    fmt_tldr = count_organism_dirs(fmt_tldr_dir, "_is_records.json") if fmt_tldr_dir.exists() else 0

    return {
        "batch": batch_name,
        "expected": expected,
        "meta_sniffles": meta_sniffles_n,
        "download": f"{dl_ok}/{dl_total}" if dl_total else "--",
        "resolve": f"{resolve_n} orgs" if resolve_n else (f"{fna_count} fna" if fna_count else "--"),
        "align": f"{align_ok}/{align_total}" if align_total else "--",
        "tldr_s": tldr_sample,
        "tldr_o": tldr_org,
        "sniffles_s": sniffles_sample,
        "sniffles_o": sniffles_org,
        "circle_s": circle_sample,
        "circle_o": circle_org,
        "format_s": fmt_sample,
        "format_o": fmt_org,
        "format_tldr": fmt_tldr,
    }


def main():
    parser = argparse.ArgumentParser(description="Check pipeline status across batches.")
    parser.add_argument("start", nargs="?", type=int, default=0, help="Start batch index")
    parser.add_argument("end", nargs="?", type=int, default=None, help="End batch index (inclusive)")
    parser.add_argument(
        "--outroot", default="/groups/rubin/projects/kuang/out/IS_cycle",
        help="Root output directory",
    )
    parser.add_argument(
        "--batchdir", default="/home/kuangh/tools/Cycle/data/batches",
        help="Directory containing batch metadata TSVs",
    )
    args = parser.parse_args()

    # Auto-detect end if not specified
    if args.end is None:
        batch_files = sorted(Path(args.batchdir).glob("batch_*.tsv"))
        if batch_files:
            args.end = int(batch_files[-1].stem.split("_")[1])
        else:
            args.end = args.start

    rows = []
    for i in range(args.start, args.end + 1):
        batch_name = f"batch_{i:03d}"
        batch_dir = Path(args.outroot) / batch_name
        if not batch_dir.exists() and not (Path(args.batchdir) / f"{batch_name}.tsv").exists():
            continue
        rows.append(check_batch(batch_dir, batch_name, args.batchdir))

    if not rows:
        print("No batches found.")
        sys.exit(1)

    # Print header
    print(f"{'Batch':<12} {'Expect':>6} {'Download':>10} {'Resolve':>12} {'Align':>10}"
          f"  {'Sniffles':>14} {'Tldr':>12} {'Circle':>14} {'Format':>14}")
    print(f"{'':12} {'':>6} {'':>10} {'':>12} {'':>10}"
          f"  {'sample/org':>14} {'sample/org':>12} {'sample/org':>14} {'sample/org':>14}")
    print("-" * 120)

    # Totals
    tot = {k: 0 for k in ["expected", "sniffles_s", "sniffles_o", "tldr_s", "tldr_o",
                           "circle_s", "circle_o", "format_s", "format_o"]}

    for r in rows:
        sniffles_str = f"{r['sniffles_s']}/{r['sniffles_o']}" if (r['sniffles_s'] or r['sniffles_o']) else "--"
        tldr_str = f"{r['tldr_s']}/{r['tldr_o']}" if (r['tldr_s'] or r['tldr_o']) else "--"
        circle_str = f"{r['circle_s']}/{r['circle_o']}" if (r['circle_s'] or r['circle_o']) else "--"
        format_str = f"{r['format_s']}/{r['format_o']}" if (r['format_s'] or r['format_o']) else "--"
        if r['format_tldr']:
            format_str += f"+{r['format_tldr']}t"

        print(f"{r['batch']:<12} {r['expected']:>6} {r['download']:>10} {r['resolve']:>12} {r['align']:>10}"
              f"  {sniffles_str:>14} {tldr_str:>12} {circle_str:>14} {format_str:>14}")

        tot["expected"] += r["expected"]
        for k in ["sniffles_s", "sniffles_o", "tldr_s", "tldr_o", "circle_s", "circle_o", "format_s", "format_o"]:
            tot[k] += r[k]

    print("-" * 120)
    print(f"{'TOTAL':<12} {tot['expected']:>6} {'':>10} {'':>12} {'':>10}"
          f"  {tot['sniffles_s']:>7}/{tot['sniffles_o']:<6}"
          f" {tot['tldr_s']:>5}/{tot['tldr_o']:<6}"
          f" {tot['circle_s']:>7}/{tot['circle_o']:<6}"
          f" {tot['format_s']:>7}/{tot['format_o']:<6}")

    print()
    print("Legend: sample/org = sample-keyed count / organism-keyed count (old pipeline)")
    print("        +Nt = N organism results from tldr-based formatter")


if __name__ == "__main__":
    main()
