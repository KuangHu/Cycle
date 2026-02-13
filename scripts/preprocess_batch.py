#!/usr/bin/env python3
"""Preprocess a batch of SRA runs into tldr-ready inputs.

Steps:
  resolve  — download one reference genome per unique organism
  index    — build minimap2 .mmi indices
  align    — minimap2 + samtools sort → sorted BAM + .bai
  is_ref   — build tldr-formatted IS element reference FASTA

Example:
    python scripts/preprocess_batch.py \\
        --metadata data/test_batch_100.tsv \\
        --fastq-dir data/sra_downloads/fastq \\
        --threads 8 \\
        --steps resolve index align is_ref
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from cycle.preprocess import ReferenceGenomeResolver, Aligner, ISReferenceBuilder, TldrRunner
from cycle.preprocess.config import (
    DEFAULT_ALIGNMENT_DIR,
    DEFAULT_FASTQ_DIR,
    DEFAULT_IS_REFERENCE_DIR,
    DEFAULT_MINIMAP2_PRESET,
    DEFAULT_REFERENCE_DIR,
    DEFAULT_SORT_MEMORY,
    DEFAULT_THREADS,
    DEFAULT_TLDR_OUTPUT_DIR,
)

ALL_STEPS = ["resolve", "index", "align", "is_ref", "tldr"]

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess SRA runs into tldr-ready inputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--metadata", required=True,
        help="Path to metadata TSV (must have 'srr_accession' and 'organism' columns).",
    )
    parser.add_argument(
        "--fastq-dir", default=DEFAULT_FASTQ_DIR,
        help=f"Directory containing FASTQ files. Default: {DEFAULT_FASTQ_DIR}",
    )
    parser.add_argument(
        "--ref-dir", default=DEFAULT_REFERENCE_DIR,
        help=f"Output directory for reference genomes. Default: {DEFAULT_REFERENCE_DIR}",
    )
    parser.add_argument(
        "--align-dir", default=DEFAULT_ALIGNMENT_DIR,
        help=f"Output directory for BAM files. Default: {DEFAULT_ALIGNMENT_DIR}",
    )
    parser.add_argument(
        "--is-dir", default=DEFAULT_IS_REFERENCE_DIR,
        help=f"Output directory for IS reference. Default: {DEFAULT_IS_REFERENCE_DIR}",
    )
    parser.add_argument(
        "--override-tsv", default=None,
        help="TSV with 'organism' and 'accession' columns for manual reference overrides.",
    )
    parser.add_argument(
        "--threads", type=int, default=DEFAULT_THREADS,
        help=f"Threads for minimap2/samtools. Default: {DEFAULT_THREADS}",
    )
    parser.add_argument(
        "--preset", default=DEFAULT_MINIMAP2_PRESET,
        help=f"minimap2 preset. Default: {DEFAULT_MINIMAP2_PRESET}",
    )
    parser.add_argument(
        "--sort-memory", default=DEFAULT_SORT_MEMORY,
        help=f"Memory per samtools sort thread. Default: {DEFAULT_SORT_MEMORY}",
    )
    parser.add_argument(
        "--is-families", nargs="*", default=None,
        help="Restrict IS reference to these families (e.g. IS256 IS6 IS3).",
    )
    parser.add_argument(
        "--tldr-dir", default=DEFAULT_TLDR_OUTPUT_DIR,
        help=f"Output directory for tldr results. Default: {DEFAULT_TLDR_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--tldr-procs", type=int, default=8,
        help="Number of processes for tldr. Default: 8",
    )
    parser.add_argument(
        "--steps", nargs="+", default=ALL_STEPS, choices=ALL_STEPS,
        help=f"Pipeline steps to run. Default: all ({' '.join(ALL_STEPS)})",
    )

    return parser.parse_args()


def find_fastq(fastq_dir: Path, sample_id: str) -> Path | None:
    """Find FASTQ file for a sample, trying common naming patterns."""
    for pattern in (
        f"{sample_id}.fastq.gz",
        f"{sample_id}.fastq",
        f"{sample_id}_1.fastq.gz",
        f"{sample_id}_pass.fastq.gz",
    ):
        path = fastq_dir / pattern
        if path.exists():
            return path

    # Glob fallback
    matches = sorted(fastq_dir.glob(f"{sample_id}*fastq*"))
    if matches:
        return matches[0]

    return None


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    metadata = pd.read_csv(args.metadata, sep="\t")
    logger.info(f"Loaded {len(metadata)} samples from {args.metadata}")

    fastq_dir = Path(args.fastq_dir)
    steps = set(args.steps)

    # ── Step 1: Resolve reference genomes ─────────────────────────────
    ref_map: dict[str, dict] = {}
    if "resolve" in steps:
        logger.info("=" * 60)
        logger.info("STEP: resolve — downloading reference genomes")
        logger.info("=" * 60)

        resolver = ReferenceGenomeResolver(
            output_dir=args.ref_dir,
            override_tsv=Path(args.override_tsv) if args.override_tsv else None,
        )
        ref_map = resolver.resolve_all(metadata)

    # ── Step 2: Build minimap2 indices ────────────────────────────────
    if "index" in steps:
        logger.info("=" * 60)
        logger.info("STEP: index — building minimap2 indices")
        logger.info("=" * 60)

        aligner = Aligner(
            output_dir=args.align_dir,
            preset=args.preset,
            threads=args.threads,
            sort_memory=args.sort_memory,
        )

        # If resolve was skipped, discover existing references
        if not ref_map:
            ref_map = _discover_existing_refs(Path(args.ref_dir))

        for org, info in ref_map.items():
            if info:
                aligner.index(info["fasta"])

    # ── Step 3: Align ─────────────────────────────────────────────────
    if "align" in steps:
        logger.info("=" * 60)
        logger.info("STEP: align — minimap2 + samtools sort")
        logger.info("=" * 60)

        aligner = Aligner(
            output_dir=args.align_dir,
            preset=args.preset,
            threads=args.threads,
            sort_memory=args.sort_memory,
        )

        if not ref_map:
            ref_map = _discover_existing_refs(Path(args.ref_dir))

        sample_map = []
        skipped = 0
        for _, row in metadata.iterrows():
            sid = row["srr_accession"]
            org = row.get("organism", "")

            ref_info = ref_map.get(org)
            if not ref_info:
                logger.warning(f"No reference for {sid} ({org}), skipping")
                skipped += 1
                continue

            fq = find_fastq(fastq_dir, sid)
            if not fq:
                logger.warning(f"No FASTQ found for {sid}, skipping")
                skipped += 1
                continue

            sample_map.append({
                "sample_id": sid,
                "fastq": fq,
                "reference_fasta": ref_info["fasta"],
            })

        if skipped:
            logger.info(f"Skipped {skipped} samples (no ref or FASTQ)")

        if sample_map:
            results = aligner.align_batch(sample_map)

            # Save alignment status
            status_path = Path(args.align_dir) / "alignment_status.tsv"
            rows = []
            for r in results:
                rows.append({
                    "sample_id": r["sample_id"],
                    "fastq": str(r["fastq"]),
                    "reference": str(r["reference_fasta"]),
                    "bam": str(r["bam"]) if r["bam"] else "",
                    "status": "ok" if r["bam"] else "failed",
                })
            pd.DataFrame(rows).to_csv(status_path, sep="\t", index=False)
            logger.info(f"Alignment status saved to {status_path}")

    # ── Step 4: IS reference ──────────────────────────────────────────
    if "is_ref" in steps:
        logger.info("=" * 60)
        logger.info("STEP: is_ref — building IS element reference")
        logger.info("=" * 60)

        builder = ISReferenceBuilder(output_dir=args.is_dir)
        builder.build(families=args.is_families)

    # ── Step 5: tldr ──────────────────────────────────────────────────
    if "tldr" in steps:
        logger.info("=" * 60)
        logger.info("STEP: tldr — detecting transposon insertions")
        logger.info("=" * 60)

        if not ref_map:
            ref_map = _discover_existing_refs(Path(args.ref_dir))

        is_ref_path = Path(args.is_dir) / "is_reference.fa"
        if not is_ref_path.exists():
            logger.error(f"IS reference not found: {is_ref_path}")
            logger.error("Run the is_ref step first.")
            sys.exit(1)

        runner = TldrRunner(
            output_dir=args.tldr_dir,
            alignment_dir=args.align_dir,
        )
        tldr_results = runner.run_batch(
            metadata=metadata,
            ref_map=ref_map,
            is_ref=is_ref_path,
            procs=args.tldr_procs,
        )

        ok = sum(1 for v in tldr_results.values() if v)
        logger.info(f"tldr: {ok}/{len(tldr_results)} organism groups produced results")

    logger.info("Done.")


def _discover_existing_refs(ref_dir: Path) -> dict[str, dict]:
    """Build ref_map from already-downloaded genomes on disk.

    This is used when the 'resolve' step was skipped. We cannot recover
    the organism→accession mapping, so we key by accession instead.
    The align step matches by organism name from metadata, so this
    helper scans existing accession dirs and returns a best-effort map.
    """
    ref_map: dict[str, dict] = {}
    if not ref_dir.exists():
        return ref_map

    for acc_dir in sorted(ref_dir.iterdir()):
        if not acc_dir.is_dir():
            continue
        fasta_files = list(acc_dir.glob("*_genomic.fna"))
        if fasta_files:
            fasta = fasta_files[0]
            fai = fasta.with_suffix(".fna.fai")
            ref_map[acc_dir.name] = {
                "accession": acc_dir.name,
                "fasta": fasta,
                "fai": fai if fai.exists() else None,
            }
    logger.info(f"Discovered {len(ref_map)} existing references in {ref_dir}")
    return ref_map


if __name__ == "__main__":
    main()
