#!/usr/bin/env python3
"""End-to-end pipeline: metadata TSV -> circle detection.

Chains up to 8 stages:
  1. download  — fetch FASTQs via kingfisher
  2. resolve   — download one reference genome per organism via NCBI datasets
  3. index     — build minimap2 .mmi indices
  4. align     — minimap2 + samtools sort -> sorted BAM + .bai
  5. is_ref    — build IS element reference FASTA (for tldr only)
  6. tldr      — run tldr per organism group (IS detection via reference library)
  6alt. sniffles — run Sniffles2 per organism (SV-based insertion detection, faster)
  7. circle    — detect IS circular intermediates via concatemer bait

Example (tldr-based):
    python scripts/run_pipeline.py --metadata batch.tsv --steps tldr circle --threads 44

Example (Sniffles2-based):
    python scripts/run_pipeline.py --metadata batch.tsv --steps sniffles circle --threads 44
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from cycle.download_manager import SRADownloader
from cycle.preprocess import ReferenceGenomeResolver, Aligner, ISReferenceBuilder, TldrRunner, SnifflesRunner
from cycle.preprocess.config import (
    DEFAULT_ALIGNMENT_DIR,
    DEFAULT_FASTQ_DIR,
    DEFAULT_IS_REFERENCE_DIR,
    DEFAULT_MINIMAP2_PRESET,
    DEFAULT_REFERENCE_DIR,
    DEFAULT_SORT_MEMORY,
    DEFAULT_THREADS,
    DEFAULT_TLDR_OUTPUT_DIR,
    DEFAULT_SNIFFLES_OUTPUT_DIR,
)
from cycle.circle_detect import CircleFinder
from cycle.circle_detect.config import (
    DEFAULT_CIRCLE_OUTPUT_DIR,
    DEFAULT_FLANK_LENGTH,
    DEFAULT_MIN_JUNCTION_OVERLAP,
)
from cycle.utils import find_fastq, slugify

ALL_STEPS = ["download", "resolve", "index", "align", "is_ref", "tldr", "sniffles", "circle"]

logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="End-to-end pipeline: metadata TSV -> circle detection.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--metadata", required=True,
        help="Path to metadata TSV (must have 'srr_accession' and 'organism' columns).",
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Root output directory.  When set, all subdirectories (fastqs, "
             "reference_genomes, alignments, …) are placed under this path "
             "and the individual --*-dir flags are ignored.",
    )
    parser.add_argument(
        "--fastq-dir", default=DEFAULT_FASTQ_DIR,
        help=f"Directory for FASTQ downloads. Default: {DEFAULT_FASTQ_DIR}",
    )
    parser.add_argument(
        "--ref-dir", default=DEFAULT_REFERENCE_DIR,
        help=f"Directory for reference genomes. Default: {DEFAULT_REFERENCE_DIR}",
    )
    parser.add_argument(
        "--align-dir", default=DEFAULT_ALIGNMENT_DIR,
        help=f"Directory for BAM files. Default: {DEFAULT_ALIGNMENT_DIR}",
    )
    parser.add_argument(
        "--is-dir", default=DEFAULT_IS_REFERENCE_DIR,
        help=f"Directory for IS reference. Default: {DEFAULT_IS_REFERENCE_DIR}",
    )
    parser.add_argument(
        "--tldr-dir", default=DEFAULT_TLDR_OUTPUT_DIR,
        help=f"Directory for tldr output. Default: {DEFAULT_TLDR_OUTPUT_DIR}",
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
        "--tldr-procs", type=int, default=8,
        help="Number of processes per tldr invocation. Default: 8",
    )
    parser.add_argument(
        "--tldr-parallel", type=int, default=1,
        help="Number of organisms to run tldr on in parallel. Default: 1",
    )
    parser.add_argument(
        "--sniffles-dir", default=None,
        help=f"Directory for Sniffles2 output. Default: derived from --outdir or {DEFAULT_SNIFFLES_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--sniffles-parallel", type=int, default=1,
        help="Number of organisms to run Sniffles2 on in parallel. Default: 1",
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
        "--steps", nargs="+", default=ALL_STEPS, choices=ALL_STEPS,
        help=f"Pipeline steps to run. Default: all ({' '.join(ALL_STEPS)})",
    )
    parser.add_argument(
        "--download-limit", type=int, default=0,
        help="Max runs to download (0 = all). Default: 0",
    )
    parser.add_argument(
        "--circle-dir", default=DEFAULT_CIRCLE_OUTPUT_DIR,
        help=f"Directory for circle detection output. Default: {DEFAULT_CIRCLE_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--min-overlap", type=int, default=DEFAULT_MIN_JUNCTION_OVERLAP,
        help=f"Min bp overlap on each side of junction. Default: {DEFAULT_MIN_JUNCTION_OVERLAP}",
    )
    parser.add_argument(
        "--flank-length", type=int, default=DEFAULT_FLANK_LENGTH,
        help=f"Genomic flank length (bp) for insertion-context baits. Default: {DEFAULT_FLANK_LENGTH}",
    )

    args = parser.parse_args()

    # When --outdir is set, derive all subdirectories from it
    if args.outdir:
        out = Path(args.outdir)
        args.fastq_dir = str(out / "sra_downloads")
        args.ref_dir = str(out / "reference_genomes")
        args.align_dir = str(out / "alignments")
        args.is_dir = str(out / "is_reference")
        args.tldr_dir = str(out / "tldr_output")
        args.sniffles_dir = str(out / "sniffles_output")
        args.circle_dir = str(out / "circle_output")

    # Set default for sniffles_dir if not provided
    if args.sniffles_dir is None:
        args.sniffles_dir = DEFAULT_SNIFFLES_OUTPUT_DIR

    return args


def discover_existing_refs(ref_dir: Path) -> dict[str, dict]:
    """Build ref_map from already-downloaded genomes on disk (accession-keyed)."""
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


def discover_organism_refs(
    align_dir: Path, metadata, ref_dir: Path,
    organism_col: str = "organism", accession_col: str = "srr_accession",
) -> dict[str, dict]:
    """Build organism-keyed ref_map from alignment_status.tsv on disk.

    Uses the alignment status file to map sample → reference FASTA, then
    joins with metadata to get organism → reference FASTA.
    """
    ref_map: dict[str, dict] = {}
    status_file = align_dir / "alignment_status.tsv"
    if not status_file.exists():
        # Fall back to accession-keyed map (won't match organism lookups,
        # but callers handle missing keys gracefully)
        return discover_existing_refs(ref_dir)

    status = pd.read_csv(status_file, sep="\t")
    # Build sample_id -> reference path
    sample_ref = {}
    for _, row in status.iterrows():
        ref_path = row.get("reference", "")
        if ref_path and row.get("status") == "ok":
            sample_ref[row["sample_id"]] = Path(ref_path)

    # Map organism -> reference FASTA (use first sample's reference per organism)
    for _, row in metadata.iterrows():
        org = row.get(organism_col, "")
        sid = row.get(accession_col, "")
        if org in ref_map or not org:
            continue
        ref_path = sample_ref.get(sid)
        if ref_path and ref_path.exists():
            fai = ref_path.with_suffix(".fna.fai")
            ref_map[org] = {
                "accession": ref_path.parent.name,
                "fasta": ref_path,
                "fai": fai if fai.exists() else None,
            }

    logger.info(f"Discovered {len(ref_map)} organism-keyed references via alignment status")
    return ref_map


def discover_tldr_tables(
    tldr_dir: Path, metadata, organism_col: str = "organism",
) -> dict[str, Path | None]:
    """Scan tldr_dir for existing .table.txt files, keyed by organism name."""
    results: dict[str, Path | None] = {}
    if not tldr_dir.exists():
        return results

    organisms = metadata[organism_col].unique()
    for organism in organisms:
        slug = slugify(organism)
        table = tldr_dir / slug / f"{slug}.table.txt"
        if table.exists():
            results[organism] = table
        else:
            results[organism] = None

    found = sum(1 for v in results.values() if v)
    logger.info(f"Discovered {found}/{len(results)} existing tldr tables in {tldr_dir}")
    return results


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
    tldr_results: dict[str, Path | None] = {}

    # ── Stage 1: Download FASTQs ─────────────────────────────────────
    if "download" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 1: download — fetching FASTQs via kingfisher")
        logger.info("=" * 60)

        downloader = SRADownloader(output_dir=args.fastq_dir)
        result_df = downloader.download_batch(
            metadata, limit=args.download_limit,
        )

        status_path = Path(args.fastq_dir) / "download_status.tsv"
        result_df.to_csv(status_path, sep="\t", index=False)
        logger.info(f"Download status saved to {status_path}")

        ok = (result_df["download_status"] == "ok").sum()
        logger.info(f"Downloaded {ok}/{len(result_df)} samples")

    # ── Stage 2: Resolve reference genomes ───────────────────────────
    ref_map: dict[str, dict] = {}
    if "resolve" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 2: resolve — downloading reference genomes")
        logger.info("=" * 60)

        resolver = ReferenceGenomeResolver(
            output_dir=args.ref_dir,
            override_tsv=Path(args.override_tsv) if args.override_tsv else None,
        )
        ref_map = resolver.resolve_all(metadata)

    # ── Stage 3: Build minimap2 indices ──────────────────────────────
    if "index" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 3: index — building minimap2 indices")
        logger.info("=" * 60)

        aligner = Aligner(
            output_dir=args.align_dir,
            preset=args.preset,
            threads=args.threads,
            sort_memory=args.sort_memory,
        )

        if not ref_map:
            ref_map = discover_existing_refs(Path(args.ref_dir))

        for org, info in ref_map.items():
            if info:
                aligner.index(info["fasta"])

    # ── Stage 4: Align reads ─────────────────────────────────────────
    if "align" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 4: align — minimap2 + samtools sort")
        logger.info("=" * 60)

        aligner = Aligner(
            output_dir=args.align_dir,
            preset=args.preset,
            threads=args.threads,
            sort_memory=args.sort_memory,
        )

        if not ref_map:
            ref_map = discover_existing_refs(Path(args.ref_dir))

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

    # ── Stage 5: IS reference ────────────────────────────────────────
    if "is_ref" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 5: is_ref — building IS element reference")
        logger.info("=" * 60)

        builder = ISReferenceBuilder(output_dir=args.is_dir)
        builder.build(families=args.is_families)

    # ── Stage 6: tldr ────────────────────────────────────────────────
    if "tldr" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 6: tldr — detecting transposon insertions")
        logger.info("=" * 60)

        if not ref_map:
            ref_map = discover_organism_refs(
                Path(args.align_dir), metadata, Path(args.ref_dir),
            )

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
            parallel=args.tldr_parallel,
        )

        ok = sum(1 for v in tldr_results.values() if v)
        logger.info(f"tldr: {ok}/{len(tldr_results)} organism groups produced results")

    # ── Stage 6alt: Sniffles2 ────────────────────────────────────────
    sniffles_results: dict[str, Path | None] = {}
    if "sniffles" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 6alt: sniffles — detecting insertions with Sniffles2")
        logger.info("=" * 60)

        if not ref_map:
            ref_map = discover_organism_refs(
                Path(args.align_dir), metadata, Path(args.ref_dir),
            )

        runner = SnifflesRunner(
            output_dir=args.sniffles_dir,
            alignment_dir=args.align_dir,
        )
        sniffles_results = runner.run_batch(
            metadata=metadata,
            ref_map=ref_map,
            parallel=args.sniffles_parallel,
        )

        ok = sum(1 for v in sniffles_results.values() if v)
        logger.info(f"Sniffles2: {ok}/{len(sniffles_results)} organism groups produced results")

    # ── Stage 7: Circle detection ─────────────────────────────────────
    if "circle" in steps:
        logger.info("=" * 60)
        logger.info("STAGE 7: circle — detecting IS circular intermediates")
        logger.info("=" * 60)

        # Discover insertion tables if neither tldr nor sniffles ran in this session
        insertion_results = {}
        if any(tldr_results.values()):
            insertion_results = tldr_results
            logger.info("Using tldr results for circle detection")
        elif any(sniffles_results.values()):
            insertion_results = sniffles_results
            logger.info("Using Sniffles2 results for circle detection")
        else:
            # Try to discover from disk
            tldr_results = discover_tldr_tables(
                Path(args.tldr_dir), metadata,
            )
            if any(tldr_results.values()):
                insertion_results = tldr_results
                logger.info("Using discovered tldr results for circle detection")
            else:
                # Try sniffles directory
                sniffles_results = discover_tldr_tables(
                    Path(args.sniffles_dir), metadata,
                )
                if any(sniffles_results.values()):
                    insertion_results = sniffles_results
                    logger.info("Using discovered Sniffles2 results for circle detection")

        if not any(insertion_results.values()):
            logger.error("No insertion tables found. Run tldr or sniffles step first.")
            sys.exit(1)

        if not ref_map:
            ref_map = discover_organism_refs(
                Path(args.align_dir), metadata, Path(args.ref_dir),
            )

        finder = CircleFinder(
            output_dir=args.circle_dir,
            min_overlap=args.min_overlap,
            flank_length=args.flank_length,
            threads=args.threads,
            sort_memory=args.sort_memory,
        )
        circle_results = finder.run_batch(
            tldr_results=insertion_results,
            metadata=metadata,
            fastq_dir=args.fastq_dir,
            ref_map=ref_map if ref_map else None,
        )

        ok = sum(1 for v in circle_results.values() if v)
        logger.info(
            f"Circle detection: {ok}/{len(circle_results)} organism groups processed"
        )

    logger.info("Pipeline complete.")


if __name__ == "__main__":
    main()
