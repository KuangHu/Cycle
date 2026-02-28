#!/usr/bin/env python3
"""Check IS_cycle data structure for consistency and completeness.

Validates per-batch directories and cross-batch analysis directories
against the canonical layout defined in DATA_STRUCTURE.md.

Usage:
    python scripts/check_data_structure.py /groups/rubin/projects/kuang/out/IS_cycle
    python scripts/check_data_structure.py /groups/rubin/projects/kuang/out/IS_cycle --batch 000
    python scripts/check_data_structure.py /groups/rubin/projects/kuang/out/IS_cycle --fix
"""

import argparse
import json
import os
import sys
from glob import glob
from pathlib import Path

# ── colours ──────────────────────────────────────────────────────────────────
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
CYAN = "\033[96m"
BOLD = "\033[1m"
RESET = "\033[0m"

OK = f"{GREEN}OK{RESET}"
WARN = f"{YELLOW}WARN{RESET}"
FAIL = f"{RED}MISSING{RESET}"
EXTRA = f"{YELLOW}EXTRA{RESET}"

# ── expected structure ───────────────────────────────────────────────────────

# Directories that belong inside a batch dir
EXPECTED_BATCH_DIRS = [
    "sra_downloads",
    "reference_genomes",
    "alignments",
    "sniffles_output",
    "circle_output",
    "is_formatter_output",
]

# Optional per-batch dirs (not an error if missing)
OPTIONAL_BATCH_DIRS = [
    "is_reference",
    "partial_circle_output",
]

# Files that may exist in a batch dir
EXPECTED_BATCH_FILES = [
    "metadata_for_sniffles.tsv",
]

OPTIONAL_BATCH_FILES = [
    "partial_circle_manifest.tsv",
]

# Known legacy/stale directories that should NOT be in a batch dir
LEGACY_DIRS = [
    "tldr_output",
    "is_formatter_output_tldr",
    "is_formatter_th_output",
    "circle_output_tldr",
    "circle_output_sniffles_test",
    "sniffles_output_test",
    "bakta_output",
    "novelty_output",  # should be at IS_cycle/ level, not inside batch
]

# Cross-batch directories at IS_cycle/ level
CROSS_BATCH_PREFIXES = [
    "system_clustering_batch_",
    "novelty_batch_",
    "is110_circular_batch_",
]

# Key files inside cross-batch dirs
CLUSTERING_EXPECTED = [
    "all_proteins.faa",
    "system_clusters.json",
    "system_clusters_summary.tsv",
]

NOVELTY_EXPECTED = [
    "blast_results.tsv",
    "cluster_novelty.json",
    "cluster_novelty_summary.tsv",
]


def fmt_size(nbytes):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if nbytes < 1024:
            return f"{nbytes:.1f}{unit}"
        nbytes /= 1024
    return f"{nbytes:.1f}PB"


def dir_size_fast(path):
    """Quick size estimate: sum immediate children sizes (no deep walk)."""
    total = 0
    try:
        for entry in os.scandir(path):
            try:
                if entry.is_file(follow_symlinks=False):
                    total += entry.stat(follow_symlinks=False).st_size
            except OSError:
                pass
    except OSError:
        pass
    return total


def count_items(path):
    try:
        return len(os.listdir(path))
    except OSError:
        return 0


# ── checks ───────────────────────────────────────────────────────────────────

def check_batch(batch_dir, batch_name, verbose=True):
    """Check a single batch directory. Returns (errors, warnings)."""
    errors = []
    warnings = []

    if verbose:
        print(f"\n{BOLD}{'=' * 60}{RESET}")
        print(f"{BOLD}  {batch_name}{RESET}")
        print(f"{BOLD}{'=' * 60}{RESET}")

    # Check expected directories
    if verbose:
        print(f"\n  {CYAN}Per-sample pipeline dirs:{RESET}")
    for d in EXPECTED_BATCH_DIRS:
        p = os.path.join(batch_dir, d)
        if os.path.isdir(p):
            n = count_items(p)
            if verbose:
                print(f"    {OK}  {d}/ ({n} items)")
            if n == 0:
                warnings.append(f"{d}/ is empty")
                if verbose:
                    print(f"         {WARN} directory is empty")
        else:
            errors.append(f"{d}/ missing")
            if verbose:
                print(f"    {FAIL}  {d}/")

    # Check optional directories
    if verbose:
        print(f"\n  {CYAN}Optional dirs:{RESET}")
    for d in OPTIONAL_BATCH_DIRS:
        p = os.path.join(batch_dir, d)
        if os.path.isdir(p):
            n = count_items(p)
            if verbose:
                print(f"    {OK}  {d}/ ({n} items)")
        else:
            if verbose:
                print(f"    --  {d}/ (not present)")

    # Check expected files
    if verbose:
        print(f"\n  {CYAN}Metadata files:{RESET}")
    for f in EXPECTED_BATCH_FILES:
        p = os.path.join(batch_dir, f)
        if os.path.isfile(p):
            if verbose:
                print(f"    {OK}  {f}")
        else:
            warnings.append(f"{f} missing")
            if verbose:
                print(f"    {WARN}  {f} (not found)")

    for f in OPTIONAL_BATCH_FILES:
        p = os.path.join(batch_dir, f)
        if os.path.isfile(p):
            if verbose:
                print(f"    {OK}  {f}")

    # Check for legacy/stale directories
    found_legacy = []
    for d in LEGACY_DIRS:
        p = os.path.join(batch_dir, d)
        if os.path.isdir(p):
            found_legacy.append(d)

    if found_legacy:
        if verbose:
            print(f"\n  {CYAN}Legacy/misplaced dirs (should be removed or moved):{RESET}")
        for d in found_legacy:
            n = count_items(os.path.join(batch_dir, d))
            warnings.append(f"legacy dir: {d}/")
            if verbose:
                print(f"    {EXTRA}  {d}/ ({n} items)")
                if d == "novelty_output":
                    print(f"         → should be at IS_cycle/novelty_{batch_name}/")

    # Check for unexpected directories
    known = set(EXPECTED_BATCH_DIRS + OPTIONAL_BATCH_DIRS + LEGACY_DIRS)
    try:
        actual = [e for e in os.listdir(batch_dir)
                  if os.path.isdir(os.path.join(batch_dir, e))]
    except OSError:
        actual = []

    unknown = [d for d in actual if d not in known]
    if unknown:
        if verbose:
            print(f"\n  {CYAN}Unexpected dirs:{RESET}")
        for d in unknown:
            n = count_items(os.path.join(batch_dir, d))
            warnings.append(f"unexpected dir: {d}/")
            if verbose:
                print(f"    {EXTRA}  {d}/ ({n} items)")

    # Check unexpected loose files
    known_files = set(EXPECTED_BATCH_FILES + OPTIONAL_BATCH_FILES)
    try:
        actual_files = [e for e in os.listdir(batch_dir)
                        if os.path.isfile(os.path.join(batch_dir, e))]
    except OSError:
        actual_files = []

    unknown_files = [f for f in actual_files if f not in known_files]
    if unknown_files:
        if verbose:
            print(f"\n  {CYAN}Unexpected files:{RESET}")
        for f in unknown_files:
            warnings.append(f"unexpected file: {f}")
            if verbose:
                print(f"    {EXTRA}  {f}")

    # Cross-check sample counts
    if verbose:
        print(f"\n  {CYAN}Sample counts:{RESET}")

    counts = {}
    for d in ["sra_downloads", "is_formatter_output", "partial_circle_output"]:
        p = os.path.join(batch_dir, d)
        if os.path.isdir(p):
            # For sra_downloads, count unique sample IDs (strip extensions)
            if d == "sra_downloads":
                items = os.listdir(p)
                samples = set()
                for item in items:
                    base = item.split(".")[0]
                    # Strip _1, _2, _subreads suffixes
                    for suffix in ("_1", "_2", "_subreads"):
                        if base.endswith(suffix):
                            base = base[: -len(suffix)]
                            break
                    samples.add(base)
                counts[d] = len(samples)
            else:
                counts[d] = count_items(p)

    for d in ["sniffles_output", "circle_output"]:
        p = os.path.join(batch_dir, d)
        if os.path.isdir(p):
            counts[d] = count_items(p)

    if verbose:
        for d, n in counts.items():
            label = "(organisms)" if d in ("sniffles_output", "circle_output") else "(samples)"
            print(f"    {d}: {n} {label}")

    # Check formatter has guide JSONs
    fmt_dir = os.path.join(batch_dir, "is_formatter_output")
    if os.path.isdir(fmt_dir):
        sample_dirs = [
            e for e in os.listdir(fmt_dir)
            if os.path.isdir(os.path.join(fmt_dir, e))
        ]
        n_guide = 0
        n_missing_guide = 0
        for sd in sample_dirs:
            guide_files = glob(os.path.join(fmt_dir, sd, "*_is_records_guide.json"))
            if guide_files:
                n_guide += 1
            else:
                n_missing_guide += 1

        if verbose:
            print(f"\n  {CYAN}Formatter completeness:{RESET}")
            print(f"    guide JSONs: {n_guide}/{len(sample_dirs)} samples")
            if n_missing_guide:
                print(f"    {WARN} {n_missing_guide} samples missing *_is_records_guide.json")
                warnings.append(
                    f"{n_missing_guide} samples missing guide JSON"
                )

    return errors, warnings


def check_cross_batch(root, verbose=True):
    """Check cross-batch analysis directories."""
    errors = []
    warnings = []

    if verbose:
        print(f"\n{BOLD}{'=' * 60}{RESET}")
        print(f"{BOLD}  Cross-batch analysis{RESET}")
        print(f"{BOLD}{'=' * 60}{RESET}")

    try:
        entries = sorted(os.listdir(root))
    except OSError:
        return errors, warnings

    # Find clustering dirs
    clustering_dirs = [
        e for e in entries if e.startswith("system_clustering_batch_")
    ]
    novelty_dirs = [e for e in entries if e.startswith("novelty_batch_")]
    is110_dirs = [e for e in entries if e.startswith("is110_circular_batch_")]

    if verbose:
        print(f"\n  {CYAN}System clustering:{RESET}")
    for d in clustering_dirs:
        p = os.path.join(root, d)
        if verbose:
            print(f"    {d}/")
        for f in CLUSTERING_EXPECTED:
            fp = os.path.join(p, f)
            if os.path.isfile(fp):
                if verbose:
                    print(f"      {OK}  {f}")
            else:
                warnings.append(f"{d}/{f} missing")
                if verbose:
                    print(f"      {FAIL}  {f}")

    if verbose:
        print(f"\n  {CYAN}Novelty annotation:{RESET}")
    for d in novelty_dirs:
        p = os.path.join(root, d)
        if verbose:
            print(f"    {d}/")
        for f in NOVELTY_EXPECTED:
            fp = os.path.join(p, f)
            if os.path.isfile(fp):
                if verbose:
                    print(f"      {OK}  {f}")
            else:
                warnings.append(f"{d}/{f} missing")
                if verbose:
                    print(f"      {FAIL}  {f}")

    if verbose and is110_dirs:
        print(f"\n  {CYAN}IS110 circular:{RESET}")
        for d in is110_dirs:
            n = count_items(os.path.join(root, d))
            print(f"    {OK}  {d}/ ({n} items)")

    # Check for clustering/novelty mismatch
    clust_batches = set()
    for d in clustering_dirs:
        tag = d.replace("system_clustering_", "")
        clust_batches.add(tag)

    nov_batches = set()
    for d in novelty_dirs:
        tag = d.replace("novelty_", "")
        nov_batches.add(tag)

    missing_novelty = clust_batches - nov_batches
    if missing_novelty and verbose:
        print(f"\n  {CYAN}Clustering without novelty:{RESET}")
        for tag in sorted(missing_novelty):
            print(f"    {WARN}  {tag} has clustering but no novelty")
            warnings.append(f"{tag}: clustering exists but novelty missing")

    # Check for novelty dirs still inside batch dirs (misplaced)
    batch_dirs = sorted(
        e for e in entries if e.startswith("batch_") and os.path.isdir(os.path.join(root, e))
    )
    misplaced = []
    for bd in batch_dirs:
        p = os.path.join(root, bd, "novelty_output")
        if os.path.isdir(p):
            misplaced.append(bd)

    if misplaced and verbose:
        print(f"\n  {CYAN}Misplaced novelty_output/ (should be at root level):{RESET}")
        for bd in misplaced:
            print(f"    {EXTRA}  {bd}/novelty_output/  →  novelty_{bd}/")
            warnings.append(f"{bd}/novelty_output should be moved to novelty_{bd}")

    return errors, warnings


def fix_misplaced_novelty(root, dry_run=True):
    """Move novelty_output/ from inside batch dirs to root level."""
    try:
        entries = sorted(os.listdir(root))
    except OSError:
        return

    batch_dirs = [
        e for e in entries
        if e.startswith("batch_") and os.path.isdir(os.path.join(root, e))
    ]

    for bd in batch_dirs:
        src = os.path.join(root, bd, "novelty_output")
        dst = os.path.join(root, f"novelty_{bd}")
        if os.path.isdir(src):
            if os.path.exists(dst):
                print(f"  SKIP  {src} → {dst} (destination exists)")
                continue
            if dry_run:
                print(f"  WOULD MOVE  {src} → {dst}")
            else:
                os.rename(src, dst)
                print(f"  MOVED  {src} → {dst}")


def main():
    parser = argparse.ArgumentParser(
        description="Check IS_cycle data structure for consistency.",
    )
    parser.add_argument(
        "root",
        help="IS_cycle root directory",
    )
    parser.add_argument(
        "--batch",
        help="Check only this batch (e.g., 000, 001)",
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Move misplaced novelty_output/ dirs to root level",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="With --fix, show what would be moved without moving",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Only print errors and warnings",
    )
    args = parser.parse_args()

    root = args.root
    if not os.path.isdir(root):
        print(f"Error: {root} is not a directory", file=sys.stderr)
        sys.exit(1)

    verbose = not args.quiet

    if args.fix:
        print(f"\n{BOLD}Fixing misplaced novelty dirs...{RESET}")
        fix_misplaced_novelty(root, dry_run=args.dry_run)
        if args.dry_run:
            print("\n(dry run — use --fix without --dry-run to actually move)")
        return

    total_errors = []
    total_warnings = []

    if args.batch:
        batch_name = f"batch_{args.batch}" if not args.batch.startswith("batch_") else args.batch
        batch_dir = os.path.join(root, batch_name)
        if not os.path.isdir(batch_dir):
            print(f"Error: {batch_dir} not found", file=sys.stderr)
            sys.exit(1)
        e, w = check_batch(batch_dir, batch_name, verbose=verbose)
        total_errors.extend(e)
        total_warnings.extend(w)
    else:
        # Check all batch dirs
        try:
            entries = sorted(os.listdir(root))
        except OSError:
            print(f"Error: cannot read {root}", file=sys.stderr)
            sys.exit(1)

        batch_dirs = [
            e for e in entries
            if e.startswith("batch_") and os.path.isdir(os.path.join(root, e))
        ]

        if not batch_dirs:
            print("No batch directories found.")
            sys.exit(1)

        if verbose:
            print(f"\n{BOLD}Found {len(batch_dirs)} batch directories{RESET}")

        for bd in batch_dirs:
            e, w = check_batch(
                os.path.join(root, bd), bd, verbose=verbose,
            )
            total_errors.extend(e)
            total_warnings.extend(w)

        # Check cross-batch dirs
        e, w = check_cross_batch(root, verbose=verbose)
        total_errors.extend(e)
        total_warnings.extend(w)

    # Summary
    print(f"\n{BOLD}{'=' * 60}{RESET}")
    print(f"{BOLD}  Summary{RESET}")
    print(f"{BOLD}{'=' * 60}{RESET}")
    print(f"  Errors:   {len(total_errors)}")
    print(f"  Warnings: {len(total_warnings)}")

    if total_errors:
        print(f"\n  {RED}Errors:{RESET}")
        for e in total_errors:
            print(f"    - {e}")
    if total_warnings:
        print(f"\n  {YELLOW}Warnings:{RESET}")
        for w in total_warnings:
            print(f"    - {w}")

    print()
    sys.exit(1 if total_errors else 0)


if __name__ == "__main__":
    main()
