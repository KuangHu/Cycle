#!/usr/bin/env python3
"""Process ECOR batch: align each sample to 4 E. coli reference genomes,
then run Sniffles + circle detection + formatter + ORF annotation + guide finder.

Usage:
    python scripts/run_ecor_batch.py \
        --fastq-dir /path/to/fastq \
        --ref-dir /path/to/reference_genomes \
        --output-dir /path/to/output \
        --threads 8
"""
import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.preprocess.aligner import Aligner
from cycle.preprocess.sniffles_runner import SnifflesRunner
from cycle.circle_detect.circle_finder import CircleFinder
from cycle.is_formatter.formatter import ISFormatter
from cycle.orf_annotator.annotator import ORFAnnotator
from cycle.guide_finder.finder import GuideFinder

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

REFERENCES = {
    "K12_MG1655": "GCF_000005845.2.fna",       # Phylogroup A
    "CFT073": "GCF_000007445.1.fna",            # Phylogroup B2
    "IAI39": "GCF_000013305.1.fna",             # Phylogroup F
    "Sakai": "GCF_000008865.2.fna",             # Phylogroup E
}


def find_fastqs(fastq_dir):
    """Find all FASTQ files in directory."""
    fastq_dir = Path(fastq_dir)
    fastqs = sorted(list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz")))
    return fastqs




def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastq-dir", required=True,
                        help="Directory containing FASTQ files")
    parser.add_argument("--ref-dir", required=True,
                        help="Directory containing reference genome FASTAs")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--references", nargs="+", default=list(REFERENCES.keys()),
                        choices=list(REFERENCES.keys()),
                        help="Which references to use (default: all 4)")
    args = parser.parse_args()

    fastqs = find_fastqs(args.fastq_dir)
    if not fastqs:
        logger.error("No FASTQ files found in %s", args.fastq_dir)
        sys.exit(1)
    logger.info("Found %d FASTQ files", len(fastqs))

    ref_dir = Path(args.ref_dir)
    out_base = Path(args.output_dir)

    for ref_name in args.references:
        ref_fasta = ref_dir / REFERENCES[ref_name]
        if not ref_fasta.exists():
            logger.error("Reference not found: %s", ref_fasta)
            continue

        logger.info("=== Reference: %s (%s) ===", ref_name, ref_fasta.name)
        ref_out = out_base / ref_name

        # Step 1: Index reference
        logger.info("Step 1: Indexing reference...")
        mmi_path = ref_fasta.with_suffix(".mmi")
        if not mmi_path.exists():
            subprocess.run([
                "minimap2", "-d", str(mmi_path), str(ref_fasta),
            ], check=True, capture_output=True)
            logger.info("Built index: %s", mmi_path)

        for fastq in fastqs:
            # Extract sample ID from filename
            sample_id = fastq.name.split(".")[0]
            logger.info("--- Sample: %s vs %s ---", sample_id, ref_name)

            sample_out = ref_out / sample_id

            # Step 2: Align
            logger.info("Step 2: Aligning...")
            align_dir = sample_out / "alignments"
            aligner = Aligner(output_dir=str(align_dir), threads=args.threads)
            bam_path = aligner.align(fastq, ref_fasta, sample_id=sample_id)
            if not bam_path:
                logger.error("Alignment failed for %s", sample_id)
                continue

            # Step 3: Sniffles (VCF + table)
            logger.info("Step 3: Sniffles...")
            sniffles_dir = sample_out / "sniffles_output"
            runner = SnifflesRunner(
                output_dir=str(sniffles_dir),
                alignment_dir=str(align_dir),
            )
            table_path = runner.run_sample(
                bam=Path(bam_path),
                ref_fasta=ref_fasta,
                sample_id=sample_id,
            )
            if not table_path:
                logger.warning("No Sniffles insertions for %s", sample_id)
                continue

            # Step 4: Circle detection
            logger.info("Step 4: Circle detection...")
            circle_dir = sample_out / "circle_output"
            try:
                cfinder = CircleFinder(output_dir=str(circle_dir), threads=args.threads)
                cfinder.run_sample(
                    tldr_table=Path(table_path),
                    sample_id=sample_id,
                    fastq_path=fastq,
                )
            except Exception as e:
                logger.error("Circle detection failed for %s: %s", sample_id, e)
                continue

            # Step 5: IS Formatter (require_th_reads=False to get all IS elements)
            logger.info("Step 5: IS Formatter...")
            circle_sample_dir = circle_dir / sample_id
            fmt_dir = sample_out / "is_formatter_output"
            os.makedirs(str(fmt_dir), exist_ok=True)
            try:
                formatter = ISFormatter(
                    output_dir=str(fmt_dir),
                    require_th_reads=False,
                )
                formatter.run_sample(
                    circle_dir=circle_sample_dir,
                    sample_id=sample_id,
                )
            except Exception as e:
                logger.error("Formatter failed for %s: %s", sample_id, e)
                continue

            # Step 6: ORF annotation
            logger.info("Step 6: ORF annotation...")
            raw_json = fmt_dir / sample_id / f"{sample_id}_is_records.json"
            ann_json = fmt_dir / sample_id / f"{sample_id}_is_records_annotated.json"
            if raw_json.exists() and not ann_json.exists():
                try:
                    annotator = ORFAnnotator()
                    annotator.annotate_sample(str(raw_json))
                except Exception as e:
                    logger.error("ORF annotation failed for %s: %s", sample_id, e)

            # Step 7: Guide finder
            logger.info("Step 7: Guide finder...")
            if ann_json.exists():
                try:
                    gf = GuideFinder(min_length=9, max_mismatches=1,
                                     min_length_for_mismatch=12)
                    gf.find_guides_sample(str(ann_json))
                except Exception as e:
                    logger.error("Guide finder failed for %s: %s", sample_id, e)

            logger.info("Done: %s vs %s", sample_id, ref_name)

    logger.info("=== All done ===")


if __name__ == "__main__":
    main()
