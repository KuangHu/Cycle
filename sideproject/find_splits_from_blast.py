"""Parse BLASTn results to find IS110 elements split by insertions.

Looks for cases where one IS110 query has two HSPs on the same subject contig,
covering different parts of the query, with a gap on the subject = inserted DNA.

Input: BLASTn outfmt 6 with columns:
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
"""
import argparse
import csv
import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import List

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    slen: int


@dataclass
class SplitResult:
    is110_id: str
    contig: str
    # Part 1
    query_start_1: int
    query_end_1: int
    ref_start_1: int
    ref_end_1: int
    pident_1: float
    # Part 2
    query_start_2: int
    query_end_2: int
    ref_start_2: int
    ref_end_2: int
    pident_2: float
    # Insertion
    insertion_start: int
    insertion_end: int
    insertion_length: int
    # Coverage
    query_length: int
    query_coverage: float


def parse_blast(path: str) -> List[BlastHit]:
    hits = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            hits.append(BlastHit(
                qseqid=parts[0], sseqid=parts[1],
                pident=float(parts[2]), length=int(parts[3]),
                qstart=int(parts[6]), qend=int(parts[7]),
                sstart=int(parts[8]), send=int(parts[9]),
                evalue=float(parts[10]), bitscore=float(parts[11]),
                qlen=int(parts[12]), slen=int(parts[13]),
            ))
    return hits


def find_splits(
    hits: List[BlastHit],
    min_query_coverage: float = 0.7,
    min_part_fraction: float = 0.15,
    min_insertion_length: int = 50,
    max_insertion_length: int = 50000,
    min_pident: float = 95.0,
    max_query_overlap: int = 50,
) -> List[SplitResult]:
    """Find pairs of HSPs on the same contig that indicate a split IS110."""

    # Group hits by (query, subject)
    by_pair = defaultdict(list)
    for h in hits:
        if h.pident >= min_pident:
            by_pair[(h.qseqid, h.sseqid)].append(h)

    results = []
    for (qid, sid), pair_hits in by_pair.items():
        if len(pair_hits) < 2:
            continue

        # Sort by query start
        pair_hits.sort(key=lambda h: h.qstart)

        # Check all pairs
        for i in range(len(pair_hits)):
            for j in range(i + 1, len(pair_hits)):
                h1, h2 = pair_hits[i], pair_hits[j]
                result = _check_pair(h1, h2, min_query_coverage, min_part_fraction,
                                     min_insertion_length, max_insertion_length,
                                     max_query_overlap)
                if result:
                    results.append(result)

    # Deduplicate: keep best per (is110_id, contig)
    best = {}
    for r in results:
        key = (r.is110_id, r.contig)
        if key not in best or r.query_coverage > best[key].query_coverage:
            best[key] = r

    return sorted(best.values(), key=lambda r: (-r.query_coverage, r.is110_id))


def _check_pair(h1, h2, min_qcov, min_part, min_ins, max_ins, max_qoverlap):
    qlen = h1.qlen

    # Query coverage of each part
    cov1 = (h1.qend - h1.qstart + 1) / qlen
    cov2 = (h2.qend - h2.qstart + 1) / qlen

    if cov1 < min_part or cov2 < min_part:
        return None

    total_cov = cov1 + cov2
    if total_cov < min_qcov:
        return None

    # Check query overlap (parts shouldn't overlap much on the query)
    q_overlap = max(0, min(h1.qend, h2.qend) - max(h1.qstart, h2.qstart) + 1)
    if q_overlap > max_qoverlap:
        return None

    # Get reference coordinates (handle reverse strand)
    rs1, re1 = min(h1.sstart, h1.send), max(h1.sstart, h1.send)
    rs2, re2 = min(h2.sstart, h2.send), max(h2.sstart, h2.send)

    # Ensure part1 is before part2 on reference
    if rs1 > rs2:
        rs1, re1, rs2, re2 = rs2, re2, rs1, re1
        h1, h2 = h2, h1
        cov1, cov2 = cov2, cov1

    # Gap on reference = insertion
    gap = rs2 - re1
    if gap < min_ins or gap > max_ins:
        return None

    return SplitResult(
        is110_id=h1.qseqid,
        contig=h1.sseqid,
        query_start_1=h1.qstart, query_end_1=h1.qend,
        ref_start_1=rs1, ref_end_1=re1, pident_1=h1.pident,
        query_start_2=h2.qstart, query_end_2=h2.qend,
        ref_start_2=rs2, ref_end_2=re2, pident_2=h2.pident,
        insertion_start=re1, insertion_end=rs2,
        insertion_length=gap,
        query_length=h1.qlen,
        query_coverage=cov1 + cov2,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast-results", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-query-coverage", type=float, default=0.7)
    parser.add_argument("--min-insertion-length", type=int, default=50)
    parser.add_argument("--max-insertion-length", type=int, default=50000)
    parser.add_argument("--min-pident", type=float, default=80.0)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    logger.info("Parsing BLAST results...")
    hits = parse_blast(args.blast_results)
    logger.info("Loaded %d BLAST hits", len(hits))

    logger.info("Finding split IS110 insertions...")
    splits = find_splits(
        hits,
        min_query_coverage=args.min_query_coverage,
        min_insertion_length=args.min_insertion_length,
        max_insertion_length=args.max_insertion_length,
        min_pident=args.min_pident,
    )
    logger.info("Found %d split IS110 insertions", len(splits))

    # Export
    tsv_path = os.path.join(args.output_dir, "is110_split_hits.tsv")
    fields = [
        "is110_id", "contig",
        "query_start_1", "query_end_1", "ref_start_1", "ref_end_1", "pident_1",
        "query_start_2", "query_end_2", "ref_start_2", "ref_end_2", "pident_2",
        "insertion_start", "insertion_end", "insertion_length",
        "query_length", "query_coverage",
    ]
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for s in splits:
            writer.writerow({k: f"{v:.1f}" if isinstance(v, float) else v
                             for k, v in s.__dict__.items()})
    logger.info("Wrote %s", tsv_path)

    # Summary stats
    unique_is110 = len(set(s.is110_id for s in splits))
    unique_contigs = len(set(s.contig for s in splits))
    ins_lengths = [s.insertion_length for s in splits]
    if ins_lengths:
        logger.info("Unique IS110 with splits: %d", unique_is110)
        logger.info("Unique contigs: %d", unique_contigs)
        logger.info("Insertion lengths: min=%d, median=%d, max=%d",
                     min(ins_lengths), sorted(ins_lengths)[len(ins_lengths)//2], max(ins_lengths))


if __name__ == "__main__":
    main()
