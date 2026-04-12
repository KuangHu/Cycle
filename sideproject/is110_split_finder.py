"""
IS110 Split Finder — detect IS110 elements with insertions in GTDB genomes.

Looks for: [IS110_part1][unknown_DNA][IS110_part2]
where IS110 = part1 + part2 (the full element is split by an insertion).

Approach:
  1. Build a concatenated GTDB genome database (or use pre-built)
  2. Map IS110 consensus sequences with minimap2 (asm20 preset for divergent seqs)
  3. Find split alignments: same IS110 query maps to same contig in two pieces
     with a gap on the reference (= inserted unknown DNA)
  4. Extract the unknown DNA between the two aligned parts

Usage:
    python sideproject/is110_split_finder.py \
        --is110-fasta /path/to/is110_consensus.fna \
        --genome-dir /path/to/gtdb_genomes/ \
        --output-dir /path/to/output \
        --threads 48
"""

import argparse
import csv
import json
import logging
import os
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import pysam

logger = logging.getLogger(__name__)


@dataclass
class SplitHit:
    """A split alignment indicating an insertion within an IS110 element."""
    is110_id: str
    genome_id: str
    contig: str
    # Part 1 alignment on reference
    ref_start_1: int
    ref_end_1: int
    query_start_1: int
    query_end_1: int
    # Part 2 alignment on reference
    ref_start_2: int
    ref_end_2: int
    query_start_2: int
    query_end_2: int
    # Insertion info
    insertion_start: int  # ref coord
    insertion_end: int    # ref coord
    insertion_length: int
    # Query coverage
    query_coverage_1: float
    query_coverage_2: float
    total_query_coverage: float
    # Alignment quality
    mapq_1: int
    mapq_2: int


class IS110SplitFinder:
    """Find IS110 elements split by insertions in genome databases."""

    def __init__(
        self,
        min_query_coverage: float = 0.7,
        max_insertion_length: int = 50000,
        min_insertion_length: int = 50,
        min_part_fraction: float = 0.15,
        min_mapq: int = 20,
        threads: int = 8,
    ):
        self.min_query_coverage = min_query_coverage
        self.max_insertion_length = max_insertion_length
        self.min_insertion_length = min_insertion_length
        self.min_part_fraction = min_part_fraction
        self.min_mapq = min_mapq
        self.threads = threads

    def build_genome_list(self, genome_dir: str) -> str:
        """Write a file listing all genome FASTA paths for minimap2."""
        genome_dir = Path(genome_dir)
        fasta_files = sorted(genome_dir.glob("*.fna"))
        if not fasta_files:
            fasta_files = sorted(genome_dir.glob("**/*.fna"))
        logger.info("Found %d genome FASTA files", len(fasta_files))
        return [str(f) for f in fasta_files]

    def run_minimap2(
        self,
        is110_fasta: str,
        genome_fastas: List[str],
        output_bam: str,
        batch_size: int = 1000,
    ) -> str:
        """Map IS110 sequences against genomes with minimap2.

        For large genome sets, concatenates in batches and maps against each.
        Uses asm20 preset for divergent sequence alignment.
        """
        os.makedirs(os.path.dirname(output_bam), exist_ok=True)

        # Process in batches to avoid too many files at once
        bam_parts = []
        for i in range(0, len(genome_fastas), batch_size):
            batch = genome_fastas[i:i + batch_size]
            batch_idx = i // batch_size
            batch_bam = output_bam.replace(".bam", f".part{batch_idx}.bam")

            if os.path.exists(batch_bam):
                logger.info("Batch %d: reusing %s", batch_idx, batch_bam)
                bam_parts.append(batch_bam)
                continue

            logger.info("Batch %d: mapping against %d genomes", batch_idx, len(batch))

            # minimap2 can take multiple reference files
            cmd = [
                "minimap2", "-a",
                "-x", "asm20",
                "-t", str(self.threads),
                "--secondary=yes",
                "-N", "1000",
                "-p", "0.5",
            ] + batch + [is110_fasta]

            sort_cmd = [
                "samtools", "sort",
                "-@", str(min(self.threads, 4)),
                "-o", batch_bam,
                "-",
            ]

            mm2 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sort = subprocess.Popen(sort_cmd, stdin=mm2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            mm2.stdout.close()
            sort.communicate()
            mm2.wait()

            subprocess.run(["samtools", "index", batch_bam], check=True)
            bam_parts.append(batch_bam)

        # Merge if multiple parts
        if len(bam_parts) == 1:
            if bam_parts[0] != output_bam:
                os.rename(bam_parts[0], output_bam)
                os.rename(bam_parts[0] + ".bai", output_bam + ".bai")
        else:
            logger.info("Merging %d BAM parts...", len(bam_parts))
            subprocess.run(
                ["samtools", "merge", "-@", str(self.threads), "-f", output_bam] + bam_parts,
                check=True,
            )
            subprocess.run(["samtools", "index", output_bam], check=True)

        return output_bam

    def find_splits(self, bam_path: str, query_lengths: dict) -> List[SplitHit]:
        """Parse BAM for split alignments indicating IS110 insertions.

        Looks for supplementary alignments where:
        - Same query maps to same contig in two pieces
        - Gap on reference between the two pieces (= insertion)
        - Combined query coverage >= min_query_coverage
        """
        splits = []
        # Group alignments by query name
        query_alns = defaultdict(list)

        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            query_alns[read.query_name].append(read)
        bam.close()

        logger.info("Processing alignments for %d queries", len(query_alns))

        for qname, alns in query_alns.items():
            qlen = query_lengths.get(qname, 0)
            if qlen == 0:
                continue

            # Group by contig
            by_contig = defaultdict(list)
            for aln in alns:
                by_contig[aln.reference_name].append(aln)

            for contig, contig_alns in by_contig.items():
                if len(contig_alns) < 2:
                    continue

                # Sort by query start
                contig_alns.sort(key=lambda a: a.query_alignment_start or 0)

                # Check all pairs for split pattern
                for i in range(len(contig_alns)):
                    for j in range(i + 1, len(contig_alns)):
                        hit = self._check_split_pair(
                            qname, qlen, contig, contig_alns[i], contig_alns[j],
                        )
                        if hit:
                            splits.append(hit)

        logger.info("Found %d split hits", len(splits))
        return splits

    def _check_split_pair(
        self, qname: str, qlen: int, contig: str,
        aln1: pysam.AlignedSegment, aln2: pysam.AlignedSegment,
    ) -> Optional[SplitHit]:
        """Check if two alignments form a valid IS110 split insertion."""
        # Get query coordinates
        qs1 = aln1.query_alignment_start or 0
        qe1 = aln1.query_alignment_end or 0
        qs2 = aln2.query_alignment_start or 0
        qe2 = aln2.query_alignment_end or 0

        # Get ref coordinates
        rs1 = aln1.reference_start
        re1 = aln1.reference_end
        rs2 = aln2.reference_start
        re2 = aln2.reference_end

        # Ensure part1 is before part2 on reference
        if rs1 > rs2:
            rs1, re1, rs2, re2 = rs2, re2, rs1, re1
            qs1, qe1, qs2, qe2 = qs2, qe2, qs1, qe1
            aln1, aln2 = aln2, aln1

        # Check for gap on reference (insertion)
        gap = rs2 - re1
        if gap < self.min_insertion_length or gap > self.max_insertion_length:
            return None

        # Check query coverage
        cov1 = (qe1 - qs1) / qlen
        cov2 = (qe2 - qs2) / qlen
        total_cov = cov1 + cov2

        if total_cov < self.min_query_coverage:
            return None

        # Each part should be substantial
        if cov1 < self.min_part_fraction or cov2 < self.min_part_fraction:
            return None

        # Check mapping quality
        if aln1.mapping_quality < self.min_mapq or aln2.mapping_quality < self.min_mapq:
            return None

        # Extract genome ID from contig name
        genome_id = contig.split()[0] if " " in contig else contig

        return SplitHit(
            is110_id=qname,
            genome_id=genome_id,
            contig=contig,
            ref_start_1=rs1, ref_end_1=re1,
            query_start_1=qs1, query_end_1=qe1,
            ref_start_2=rs2, ref_end_2=re2,
            query_start_2=qs2, query_end_2=qe2,
            insertion_start=re1,
            insertion_end=rs2,
            insertion_length=gap,
            query_coverage_1=cov1,
            query_coverage_2=cov2,
            total_query_coverage=total_cov,
            mapq_1=aln1.mapping_quality,
            mapq_2=aln2.mapping_quality,
        )

    def export_results(self, splits: List[SplitHit], output_dir: str):
        """Write split hits as TSV."""
        os.makedirs(output_dir, exist_ok=True)
        tsv_path = os.path.join(output_dir, "is110_split_hits.tsv")

        fields = [
            "is110_id", "genome_id", "contig",
            "ref_start_1", "ref_end_1", "query_start_1", "query_end_1",
            "ref_start_2", "ref_end_2", "query_start_2", "query_end_2",
            "insertion_start", "insertion_end", "insertion_length",
            "query_coverage_1", "query_coverage_2", "total_query_coverage",
            "mapq_1", "mapq_2",
        ]

        with open(tsv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for hit in splits:
                writer.writerow({
                    k: f"{v:.3f}" if isinstance(v, float) else v
                    for k, v in hit.__dict__.items()
                })

        logger.info("Wrote %d split hits to %s", len(splits), tsv_path)
        return tsv_path

    def run(
        self,
        is110_fasta: str,
        genome_dir: str,
        output_dir: str,
    ) -> List[SplitHit]:
        """Full pipeline: map → find splits → export."""
        os.makedirs(output_dir, exist_ok=True)

        # Get query lengths
        query_lengths = {}
        with open(is110_fasta) as f:
            name = None
            seq_len = 0
            for line in f:
                if line.startswith(">"):
                    if name:
                        query_lengths[name] = seq_len
                    name = line[1:].strip().split()[0]
                    seq_len = 0
                else:
                    seq_len += len(line.strip())
            if name:
                query_lengths[name] = seq_len
        logger.info("Loaded %d IS110 query sequences", len(query_lengths))

        # Get genome list
        genome_fastas = self.build_genome_list(genome_dir)

        # Map
        bam_path = os.path.join(output_dir, "is110_vs_gtdb.bam")
        self.run_minimap2(is110_fasta, genome_fastas, bam_path)

        # Find splits
        splits = self.find_splits(bam_path, query_lengths)

        # Export
        self.export_results(splits, output_dir)

        return splits


def main():
    parser = argparse.ArgumentParser(
        description="Find IS110 elements split by insertions in GTDB genomes.",
    )
    parser.add_argument("--is110-fasta", required=True,
                        help="FASTA of IS110 consensus sequences")
    parser.add_argument("--genome-dir", required=True,
                        help="Directory with genome FASTA files (*.fna)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-query-coverage", type=float, default=0.7)
    parser.add_argument("--max-insertion-length", type=int, default=50000)
    parser.add_argument("--min-insertion-length", type=int, default=50)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=1000,
                        help="Genomes per minimap2 batch")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    finder = IS110SplitFinder(
        min_query_coverage=args.min_query_coverage,
        max_insertion_length=args.max_insertion_length,
        min_insertion_length=args.min_insertion_length,
        min_mapq=args.min_mapq,
        threads=args.threads,
    )
    finder.run(
        is110_fasta=args.is110_fasta,
        genome_dir=args.genome_dir,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
