"""Detect partial IS circular intermediates via split-read back-jump analysis.

A partial circle ``[S, E]`` within an IS element produces chimeric reads whose
supplementary alignments jump *backward* on the linear IS reference — from
position E back to position S.  This module detects those back-jumps, clusters
nearby breakpoints, and reports partial circle calls per IS element.
"""

import csv
import json
import logging
import math
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Optional

import pysam

from .config import (
    DEFAULT_BREAKPOINT_TOLERANCE,
    DEFAULT_MAX_CIRCLE_FRACTION,
    DEFAULT_MIN_CIRCLE_SIZE,
    DEFAULT_MIN_CONSENSUS_ENTROPY,
    DEFAULT_MIN_CONSENSUS_LENGTH,
    DEFAULT_MIN_OVERLAP_EACH_SIDE,
    DEFAULT_MIN_SUPPORTING_READS,
    DEFAULT_PARTIAL_CIRCLE_OUTPUT_DIR,
)

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# Data structures
# ------------------------------------------------------------------

@dataclass
class SplitSegment:
    """One alignment segment (primary or supplementary) of a split read."""

    read_name: str
    is_id: str
    ref_start: int      # 0-based inclusive
    ref_end: int         # 0-based exclusive
    query_start: int     # position on read (oriented by strand)
    query_end: int
    is_reverse: bool
    mapping_quality: int


@dataclass
class BackJump:
    """A detected back-jump from one split segment to the next."""

    read_name: str
    is_id: str
    circle_start: int    # S — where the read jumps TO (ref_start of later segment)
    circle_end: int      # E — where the read jumps FROM (ref_end of earlier segment)
    circle_size: int
    overlap_before: int  # bp aligned before the junction
    overlap_after: int   # bp aligned after the junction
    mapping_quality: int


@dataclass
class PartialCircleCall:
    """A clustered partial circle call with supporting reads."""

    is_id: str
    is_length: int
    circle_start: int
    circle_end: int
    circle_size: int
    circle_fraction: float
    supporting_reads: list[BackJump] = field(default_factory=list)

    @property
    def n_supporting_reads(self) -> int:
        return len(self.supporting_reads)

    @property
    def mean_mapq(self) -> float:
        if not self.supporting_reads:
            return 0.0
        return sum(r.mapping_quality for r in self.supporting_reads) / len(
            self.supporting_reads
        )


class PartialCircleDetector:
    """Detect partial IS circular intermediates via split-read analysis.

    Maps reads to single-copy IS references, identifies supplementary
    alignments that jump backward on the reference (back-jumps), clusters
    nearby breakpoints, and reports partial circle calls.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_PARTIAL_CIRCLE_OUTPUT_DIR,
        min_overlap: int = DEFAULT_MIN_OVERLAP_EACH_SIDE,
        min_circle_size: int = DEFAULT_MIN_CIRCLE_SIZE,
        max_circle_fraction: float = DEFAULT_MAX_CIRCLE_FRACTION,
        breakpoint_tolerance: int = DEFAULT_BREAKPOINT_TOLERANCE,
        min_supporting_reads: int = DEFAULT_MIN_SUPPORTING_READS,
        min_consensus_length: int = DEFAULT_MIN_CONSENSUS_LENGTH,
        min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
        threads: int = 8,
        sort_memory: str = "4G",
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_overlap = min_overlap
        self.min_circle_size = min_circle_size
        self.max_circle_fraction = max_circle_fraction
        self.breakpoint_tolerance = breakpoint_tolerance
        self.min_supporting_reads = min_supporting_reads
        self.min_consensus_length = min_consensus_length
        self.min_consensus_entropy = min_consensus_entropy
        self.threads = threads
        self.sort_memory = sort_memory

        for tool in ("minimap2", "samtools"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _seq_entropy(seq: str) -> float:
        """Shannon entropy in bits per base (0 = homopolymer, 2 = random)."""
        if not seq:
            return 0.0
        seq = seq.upper()
        n = len(seq)
        counts = Counter(seq)
        return -sum((c / n) * math.log2(c / n) for c in counts.values())

    def _build_is_reference(
        self, entries: list[dict], output_path: Path,
    ) -> tuple[Path, dict[str, int]]:
        """Write single-copy IS reference FASTA — one entry per is_id.

        Args:
            entries: List of IS record dicts with ``is_id`` and
                ``is_element.sequence`` fields.
            output_path: Where to write the FASTA.

        Returns:
            ``(fasta_path, is_lengths)`` where *is_lengths* maps is_id to
            sequence length.
        """
        if output_path.exists():
            logger.info(f"IS reference FASTA exists: {output_path}")
            is_lengths: dict[str, int] = {}
            with open(output_path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        is_id = line[1:].strip().split()[0]
                    else:
                        is_lengths[is_id] = is_lengths.get(is_id, 0) + len(
                            line.strip()
                        )
            return output_path, is_lengths

        is_lengths = {}
        written = 0
        skipped_short = 0
        skipped_low_complexity = 0

        with open(output_path, "w") as fh:
            for rec in entries:
                is_id = rec["is_id"]
                seq = (rec.get("is_element") or {}).get("sequence", "")
                if not seq:
                    continue
                if len(seq) < self.min_consensus_length:
                    skipped_short += 1
                    continue
                if self._seq_entropy(seq) < self.min_consensus_entropy:
                    skipped_low_complexity += 1
                    continue

                fh.write(f">{is_id}\n{seq}\n")
                is_lengths[is_id] = len(seq)
                written += 1

        logger.info(
            f"Wrote {written} IS sequences to {output_path.name}"
            f" (skipped {skipped_short} short, "
            f"{skipped_low_complexity} low-complexity)"
        )
        return output_path, is_lengths

    def _map_reads(
        self, fastq: Path, ref_fasta: Path, output_bam: Path,
    ) -> Optional[Path]:
        """Map reads to single-copy IS reference with minimap2.

        Uses name-sorted output (``samtools sort -n``) for efficient
        single-pass grouping of primary + supplementary segments by read.
        Supplementary alignments are enabled by default in minimap2.
        """
        fastq = Path(fastq)
        output_bam = Path(output_bam)

        if output_bam.exists():
            logger.info(f"Name-sorted BAM exists: {output_bam}")
            return output_bam

        mm2_cmd = [
            "minimap2",
            "-a",
            "-x", "map-ont",
            "-t", str(self.threads),
            str(ref_fasta),
            str(fastq),
        ]
        sort_cmd = [
            "samtools", "sort",
            "-n",  # name-sorted for split-read grouping
            "-@", str(self.threads),
            "-m", self.sort_memory,
            "-o", str(output_bam),
            "-",
        ]

        logger.info(f"Mapping {fastq.name} -> {output_bam.name} (name-sorted)")
        try:
            mm2_proc = subprocess.Popen(
                mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            sort_proc = subprocess.Popen(
                sort_cmd,
                stdin=mm2_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            mm2_proc.stdout.close()

            sort_out, sort_err = sort_proc.communicate()
            mm2_err = mm2_proc.stderr.read()
            mm2_proc.stderr.close()
            mm2_proc.wait()

            if sort_proc.returncode != 0:
                logger.error(
                    f"samtools sort failed: {sort_err.decode().strip()}"
                )
                return None

        except Exception as exc:
            logger.error(f"Mapping pipeline failed: {exc}")
            return None

        logger.info(f"  -> {output_bam}")
        return output_bam

    def _extract_split_segments(
        self, bam_path: Path, is_id_set: set[str],
    ) -> dict[str, list[SplitSegment]]:
        """Iterate name-sorted BAM and group segments by (read_name, is_id).

        Collects both primary and supplementary alignments.  Secondary
        alignments and unmapped reads are skipped.

        Returns:
            Dict mapping ``read_name`` to list of SplitSegments.
        """
        segments_by_read: dict[str, list[SplitSegment]] = {}

        bam = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary:
                continue

            ref_name = read.reference_name
            if not ref_name or ref_name not in is_id_set:
                continue

            # Get query coordinates (how much of the read is aligned)
            pairs = read.get_aligned_pairs(matches_only=False)
            query_positions = [
                qpos for qpos, rpos in pairs if qpos is not None and rpos is not None
            ]
            if not query_positions:
                continue

            seg = SplitSegment(
                read_name=read.query_name,
                is_id=ref_name,
                ref_start=read.reference_start,
                ref_end=read.reference_end,
                query_start=min(query_positions),
                query_end=max(query_positions) + 1,
                is_reverse=read.is_reverse,
                mapping_quality=read.mapping_quality,
            )

            key = read.query_name
            if key not in segments_by_read:
                segments_by_read[key] = []
            segments_by_read[key].append(seg)

        bam.close()

        # Keep only reads with 2+ segments on the same IS
        split_reads = {}
        for read_name, segs in segments_by_read.items():
            if len(segs) < 2:
                continue
            split_reads[read_name] = segs

        logger.info(
            f"Found {len(split_reads)} split reads with 2+ segments"
        )
        return split_reads

    def _detect_back_jumps(
        self,
        segments_by_read: dict[str, list[SplitSegment]],
        is_lengths: dict[str, int],
    ) -> list[BackJump]:
        """Detect back-jumps in split reads.

        For each read, group segments by (is_id, strand), sort by query
        position, and check consecutive pairs for reference coordinate
        back-jumps: ``seg_i.ref_end > seg_{i+1}.ref_start``.
        """
        back_jumps: list[BackJump] = []

        for read_name, segs in segments_by_read.items():
            # Group by (is_id, strand)
            groups: dict[tuple[str, bool], list[SplitSegment]] = {}
            for seg in segs:
                key = (seg.is_id, seg.is_reverse)
                if key not in groups:
                    groups[key] = []
                groups[key].append(seg)

            for (is_id, _strand), group in groups.items():
                if len(group) < 2:
                    continue

                is_len = is_lengths.get(is_id)
                if is_len is None:
                    continue

                # Sort by query position (read order)
                group.sort(key=lambda s: s.query_start)

                for i in range(len(group) - 1):
                    seg_a = group[i]
                    seg_b = group[i + 1]

                    # Back-jump: seg_a aligns to later ref position,
                    # seg_b aligns to earlier ref position
                    if seg_a.ref_end <= seg_b.ref_start:
                        continue  # forward jump — normal, not a back-jump

                    # Circle = [seg_b.ref_start, seg_a.ref_end]
                    circle_start = seg_b.ref_start
                    circle_end = seg_a.ref_end
                    circle_size = circle_end - circle_start

                    if circle_size < self.min_circle_size:
                        continue

                    circle_fraction = circle_size / is_len
                    if circle_fraction > self.max_circle_fraction:
                        continue

                    # Check overlap on each side of junction
                    overlap_before = seg_a.ref_end - seg_a.ref_start
                    overlap_after = seg_b.ref_end - seg_b.ref_start
                    if (overlap_before < self.min_overlap
                            or overlap_after < self.min_overlap):
                        continue

                    mapq = min(seg_a.mapping_quality, seg_b.mapping_quality)

                    back_jumps.append(BackJump(
                        read_name=read_name,
                        is_id=is_id,
                        circle_start=circle_start,
                        circle_end=circle_end,
                        circle_size=circle_size,
                        overlap_before=overlap_before,
                        overlap_after=overlap_after,
                        mapping_quality=mapq,
                    ))

        logger.info(f"Detected {len(back_jumps)} back-jumps")
        return back_jumps

    def _cluster_breakpoints(
        self,
        back_jumps: list[BackJump],
        is_lengths: dict[str, int],
    ) -> list[PartialCircleCall]:
        """Cluster back-jumps within ±tolerance and take median coordinates.

        Groups back-jumps by is_id, then uses single-linkage clustering
        on (circle_start, circle_end) within the tolerance window.
        """
        # Group by is_id
        by_is: dict[str, list[BackJump]] = {}
        for bj in back_jumps:
            if bj.is_id not in by_is:
                by_is[bj.is_id] = []
            by_is[bj.is_id].append(bj)

        calls: list[PartialCircleCall] = []
        tol = self.breakpoint_tolerance

        for is_id, jumps in by_is.items():
            is_len = is_lengths.get(is_id, 0)

            # Sort by (circle_start, circle_end)
            jumps.sort(key=lambda j: (j.circle_start, j.circle_end))

            # Single-linkage clustering
            clusters: list[list[BackJump]] = []
            for bj in jumps:
                merged = False
                for cluster in clusters:
                    # Check if this back-jump is close to any member
                    for member in cluster:
                        if (abs(bj.circle_start - member.circle_start) <= tol
                                and abs(bj.circle_end - member.circle_end) <= tol):
                            cluster.append(bj)
                            merged = True
                            break
                    if merged:
                        break
                if not merged:
                    clusters.append([bj])

            for cluster in clusters:
                if len(cluster) < self.min_supporting_reads:
                    continue

                med_start = int(median(bj.circle_start for bj in cluster))
                med_end = int(median(bj.circle_end for bj in cluster))
                circle_size = med_end - med_start

                if circle_size < self.min_circle_size:
                    continue

                circle_fraction = circle_size / is_len if is_len > 0 else 0.0
                if circle_fraction > self.max_circle_fraction:
                    continue

                calls.append(PartialCircleCall(
                    is_id=is_id,
                    is_length=is_len,
                    circle_start=med_start,
                    circle_end=med_end,
                    circle_size=circle_size,
                    circle_fraction=round(circle_fraction, 4),
                    supporting_reads=cluster,
                ))

        logger.info(
            f"Clustered into {len(calls)} partial circle calls "
            f"(min {self.min_supporting_reads} reads)"
        )
        return calls

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_sample(
        self,
        is_records_path: Path,
        sample_id: str,
        fastq_path: Path,
    ) -> Optional[Path]:
        """Run partial circle detection for one sample.

        Args:
            is_records_path: Path to ``*_is_records_guide.json`` containing
                IS element records with sequences.
            sample_id: Sample accession (e.g. SRR123456).
            fastq_path: Path to the sample's FASTQ file.

        Returns:
            Path to the sample output directory, or ``None`` on failure.
        """
        is_records_path = Path(is_records_path)
        fastq_path = Path(fastq_path)

        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 1. Load IS records and extract sequences to FASTA
        try:
            with open(is_records_path) as fh:
                records = json.load(fh)
        except Exception as exc:
            logger.error(f"Failed to load IS records from {is_records_path}: {exc}")
            return None

        if not records:
            logger.warning(f"No IS records for {sample_id}, skipping")
            return None

        ref_fasta = sample_dir / f"{sample_id}_is_ref.fa"
        ref_fasta, is_lengths = self._build_is_reference(records, ref_fasta)

        if not is_lengths:
            logger.warning(
                f"No valid IS sequences for {sample_id} after filtering, skipping"
            )
            return None

        is_id_set = set(is_lengths.keys())
        logger.info(
            f"{sample_id}: {len(is_id_set)} IS elements, mapping reads..."
        )

        # 2. Map reads with minimap2 (name-sorted BAM)
        bam_path = sample_dir / f"{sample_id}.partial.nsorted.bam"
        bam = self._map_reads(fastq_path, ref_fasta, bam_path)
        if bam is None:
            logger.warning(f"Mapping failed for {sample_id}, skipping")
            return None

        # 3. Extract split segments
        segments_by_read = self._extract_split_segments(bam, is_id_set)
        if not segments_by_read:
            logger.info(f"No split reads for {sample_id}")
            self._write_empty_outputs(sample_dir, sample_id)
            return sample_dir

        # 4. Detect back-jumps
        back_jumps = self._detect_back_jumps(segments_by_read, is_lengths)
        if not back_jumps:
            logger.info(f"No back-jumps detected for {sample_id}")
            self._write_empty_outputs(sample_dir, sample_id)
            return sample_dir

        # 5. Cluster breakpoints and filter
        calls = self._cluster_breakpoints(back_jumps, is_lengths)

        # 6. Write outputs
        self._write_reads_tsv(sample_dir, sample_id, back_jumps)
        self._write_summary_json(sample_dir, sample_id, calls)
        self._write_summary_tsv(sample_dir, sample_id, calls)

        logger.info(
            f"Partial circle detection for {sample_id}: "
            f"{len(back_jumps)} back-jumps, {len(calls)} calls"
        )
        return sample_dir

    # ------------------------------------------------------------------
    # Output writers
    # ------------------------------------------------------------------

    def _write_empty_outputs(self, sample_dir: Path, sample_id: str) -> None:
        """Write empty output files when no results are found."""
        self._write_reads_tsv(sample_dir, sample_id, [])
        self._write_summary_json(sample_dir, sample_id, [])
        self._write_summary_tsv(sample_dir, sample_id, [])

    def _write_reads_tsv(
        self, sample_dir: Path, sample_id: str, back_jumps: list[BackJump],
    ) -> None:
        """Write per-read back-jump TSV."""
        path = sample_dir / f"{sample_id}_partial_circle_reads.tsv"
        columns = [
            "read_name", "sample_id", "is_id", "circle_start", "circle_end",
            "circle_size", "overlap_before", "overlap_after", "mapping_quality",
        ]

        rows = []
        for bj in back_jumps:
            rows.append({
                "read_name": bj.read_name,
                "sample_id": sample_id,
                "is_id": bj.is_id,
                "circle_start": bj.circle_start,
                "circle_end": bj.circle_end,
                "circle_size": bj.circle_size,
                "overlap_before": bj.overlap_before,
                "overlap_after": bj.overlap_after,
                "mapping_quality": bj.mapping_quality,
            })

        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

        logger.info(f"Wrote {len(rows)} rows to {path.name}")

    def _write_summary_json(
        self, sample_dir: Path, sample_id: str, calls: list[PartialCircleCall],
    ) -> None:
        """Write per-IS partial circle summary JSON."""
        path = sample_dir / f"{sample_id}_partial_circle_summary.json"

        records = []
        for call in calls:
            records.append({
                "is_id": call.is_id,
                "sample_id": sample_id,
                "is_length": call.is_length,
                "circle_start": call.circle_start,
                "circle_end": call.circle_end,
                "circle_size": call.circle_size,
                "circle_fraction": call.circle_fraction,
                "n_supporting_reads": call.n_supporting_reads,
                "mean_mapq": round(call.mean_mapq, 1),
                "supporting_read_names": [
                    r.read_name for r in call.supporting_reads
                ],
            })

        with open(path, "w") as fh:
            json.dump(records, fh, indent=2)

        logger.info(f"Wrote {len(records)} calls to {path.name}")

    def _write_summary_tsv(
        self, sample_dir: Path, sample_id: str, calls: list[PartialCircleCall],
    ) -> None:
        """Write flat one-row-per-call summary TSV."""
        path = sample_dir / f"{sample_id}_partial_circle_summary.tsv"
        columns = [
            "is_id", "sample_id", "is_length", "circle_start", "circle_end",
            "circle_size", "circle_fraction", "n_supporting_reads", "mean_mapq",
        ]

        rows = []
        for call in calls:
            rows.append({
                "is_id": call.is_id,
                "sample_id": sample_id,
                "is_length": call.is_length,
                "circle_start": call.circle_start,
                "circle_end": call.circle_end,
                "circle_size": call.circle_size,
                "circle_fraction": call.circle_fraction,
                "n_supporting_reads": call.n_supporting_reads,
                "mean_mapq": round(call.mean_mapq, 1),
            })

        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

        logger.info(f"Wrote {len(rows)} rows to {path.name}")
