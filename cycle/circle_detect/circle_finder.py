"""Detect IS circular intermediates and insertion boundaries via bait mapping.

For each IS consensus from tldr, a concatemer bait ``[IS][IS]`` is built.
Reads from circular intermediates map across the tail-head junction at the
center seam.  Genome-head and tail-genome junctions are detected via
soft-clip analysis on reads aligned near IS copy boundaries within the
same bait — no reference genome flanking sequence is needed.
"""

import csv
import logging
import math
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pysam

from ..utils import find_fastq
from .config import (
    DEFAULT_BOUNDARY_TOLERANCE,
    DEFAULT_CIRCLE_OUTPUT_DIR,
    DEFAULT_MIN_CONSENSUS_ENTROPY,
    DEFAULT_MIN_CONSENSUS_LENGTH,
    DEFAULT_MIN_JUNCTION_OVERLAP,
)

logger = logging.getLogger(__name__)


@dataclass
class ISEntry:
    """A single IS element insertion from a tldr table."""

    uuid: str
    chrom: str
    start: int
    end: int
    family: str
    subfamily: str
    consensus: str


class CircleFinder:
    """Detect IS circular intermediates and insertion boundaries via bait mapping.

    For each sample, concatemer bait ``[IS][IS]`` sequences are built for all
    IS elements from tldr/sniffles.  Tail-head junctions are detected by reads
    spanning the center seam.  Genome-head and tail-genome junctions are
    detected via soft-clip analysis on reads aligned near IS copy boundaries.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_CIRCLE_OUTPUT_DIR,
        min_overlap: int = DEFAULT_MIN_JUNCTION_OVERLAP,
        min_consensus_length: int = DEFAULT_MIN_CONSENSUS_LENGTH,
        boundary_tolerance: int = DEFAULT_BOUNDARY_TOLERANCE,
        min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
        threads: int = 8,
        sort_memory: str = "4G",
        flank_length: int | None = None,  # deprecated, ignored
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_overlap = min_overlap
        self.min_consensus_length = min_consensus_length
        self.boundary_tolerance = boundary_tolerance
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

    def _parse_tldr_table(self, table_path: Path) -> list[ISEntry]:
        """Parse a tldr .table.txt and return IS entries with valid consensus.

        Filters out entries whose consensus is empty, 'NA', shorter than
        ``min_consensus_length``, or low-complexity (entropy below
        ``min_consensus_entropy``).
        """
        table_path = Path(table_path)
        entries: list[ISEntry] = []
        skipped_low_complexity = 0

        with open(table_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                consensus = row.get("Consensus", "").strip()
                if not consensus or consensus == "NA":
                    continue
                if len(consensus) < self.min_consensus_length:
                    continue
                if self._seq_entropy(consensus) < self.min_consensus_entropy:
                    skipped_low_complexity += 1
                    continue

                entries.append(ISEntry(
                    uuid=row["UUID"],
                    chrom=row["Chrom"],
                    start=int(row["Start"]),
                    end=int(row["End"]),
                    family=row.get("Family", "NA"),
                    subfamily=row.get("Subfamily", "NA"),
                    consensus=consensus,
                ))

        logger.info(
            f"Parsed {len(entries)} IS entries with consensus >= "
            f"{self.min_consensus_length} bp from {table_path.name}"
            f" (skipped {skipped_low_complexity} low-complexity)"
        )
        return entries

    def _build_bait_fasta(
        self,
        entries: list[ISEntry],
        output_path: Path,
    ) -> tuple[Path, dict[str, list[dict]]]:
        """Write concatemer bait FASTA ``[IS][IS]`` for each IS element.

        Header format::

            >{uuid}__th__j{N}

        Returns:
            ``(fasta_path, junction_map)`` where *junction_map* maps
            ref_name to a list of ``{type, position, uuid}`` dicts.
        """
        output_path = Path(output_path)
        if output_path.exists():
            logger.info(f"Bait FASTA exists: {output_path}")
            return output_path, self._parse_bait_headers(output_path)

        junction_map: dict[str, list[dict]] = {}

        with open(output_path, "w") as fh:
            for entry in entries:
                n = len(entry.consensus)

                th_name = f"{entry.uuid}__th__j{n}"
                fh.write(f">{th_name}\n")
                fh.write(f"{entry.consensus}{entry.consensus}\n")
                junction_map[th_name] = [
                    {"type": "tail_head", "position": n, "uuid": entry.uuid},
                ]

        logger.info(
            f"Wrote bait FASTA with {len(junction_map)} tail-head "
            f"sequences to {output_path}"
        )
        return output_path, junction_map

    def _parse_bait_headers(self, fasta_path: Path) -> dict[str, list[dict]]:
        """Reconstruct junction_map from existing bait FASTA headers.

        Handles current ``__th__`` format as well as the legacy ``__len{N}``
        format for backward compatibility.  Old ``__ic__`` headers are ignored.
        """
        junction_map: dict[str, list[dict]] = {}

        with open(fasta_path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                ref_name = line[1:].strip().split()[0]

                if "__th__j" in ref_name:
                    parts = ref_name.split("__th__j")
                    uuid = parts[0]
                    junction_pos = int(parts[1])
                    junction_map[ref_name] = [
                        {"type": "tail_head", "position": junction_pos, "uuid": uuid},
                    ]

                elif "__len" in ref_name:
                    parts = ref_name.rsplit("__len", 1)
                    uuid = parts[0]
                    junction_pos = int(parts[1])
                    junction_map[ref_name] = [
                        {"type": "tail_head", "position": junction_pos, "uuid": uuid},
                    ]

        return junction_map

    def _map_reads(
        self, fastq: Path, bait_fasta: Path, output_bam: Path,
    ) -> Optional[Path]:
        """Map reads to the bait reference with minimap2.

        Uses ``--secondary=no`` so each read maps to its best bait target only.
        Pipes minimap2 directly into samtools sort (no intermediate SAM).
        """
        fastq = Path(fastq)
        output_bam = Path(output_bam)

        if output_bam.exists() and Path(str(output_bam) + ".bai").exists():
            logger.info(f"BAM exists: {output_bam}")
            return output_bam

        mm2_cmd = [
            "minimap2",
            "-a",
            "--secondary=no",
            "-x", "map-ont",
            "-t", str(self.threads),
            str(bait_fasta),
            str(fastq),
        ]
        sort_cmd = [
            "samtools", "sort",
            "-@", str(self.threads),
            "-m", self.sort_memory,
            "-o", str(output_bam),
            "-",
        ]

        logger.info(f"Mapping {fastq.name} -> {output_bam.name}")
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

        # Index the BAM
        ret = subprocess.run(
            ["samtools", "index", str(output_bam)],
            capture_output=True, text=True,
        )
        if ret.returncode != 0:
            logger.error(f"samtools index failed: {ret.stderr.strip()}")
            return None

        logger.info(f"  -> {output_bam} + .bai")
        return output_bam

    def _find_junction_reads(
        self,
        bam_path: Path,
        entries_by_uuid: dict[str, ISEntry],
        junction_map: dict[str, list[dict]],
        sample_id: str,
    ) -> tuple[list[dict], dict[str, dict]]:
        """Scan BAM for reads spanning bait junctions.

        **Tail-head**: read alignment spans the junction at position N
        (center of ``[IS][IS]`` bait) by at least ``min_overlap`` on each side.

        **Genome-head**: read has a large left soft-clip (genomic) and its
        alignment starts near an IS copy boundary (0 or N on the TH bait).

        **Tail-genome**: read has a large right soft-clip (genomic) and its
        alignment ends near an IS copy boundary (N or 2N on the TH bait).

        Returns:
            ``(junction_reads, summary_by_uuid)`` where *junction_reads* is a
            list of per-read dicts, and *summary_by_uuid* maps uuid to counts
            split by junction type.
        """
        junction_reads: list[dict] = []
        total_mapped: dict[str, int] = {}
        counts: dict[str, dict[str, int]] = {}  # uuid -> {type -> count}
        # uuid -> (mapq, length, sequence) for best tail-head example read
        th_example: dict[str, tuple[int, int, str]] = {}

        tol = self.boundary_tolerance
        min_overlap = self.min_overlap

        bam = pysam.AlignmentFile(str(bam_path), "rb")
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            ref_name = read.reference_name
            if not ref_name or ref_name not in junction_map:
                continue

            junctions = junction_map[ref_name]
            uuid = junctions[0]["uuid"]

            entry = entries_by_uuid.get(uuid)
            if entry is None:
                continue

            total_mapped[uuid] = total_mapped.get(uuid, 0) + 1

            ref_start = read.reference_start   # 0-based
            ref_end = read.reference_end       # 0-based, exclusive
            aligned_len = ref_end - ref_start

            # Soft-clip lengths from CIGAR
            cigar = read.cigartuples
            left_clip = cigar[0][1] if cigar and cigar[0][0] in (4, 5) else 0
            right_clip = cigar[-1][1] if cigar and cigar[-1][0] in (4, 5) else 0

            # The junction_map for a TH bait has one entry with position = N
            junction_pos = junctions[0]["position"]
            n = junction_pos  # IS consensus length

            detected_types: list[str] = []

            # Tail-head: alignment spans the junction at N
            if (ref_start <= junction_pos - min_overlap
                    and ref_end >= junction_pos + min_overlap):
                detected_types.append("tail_head")

            # Genome-head: left soft-clip is genomic, alignment starts
            # near an IS copy boundary (position 0 or N)
            if (left_clip >= min_overlap
                    and aligned_len >= min_overlap
                    and (ref_start < tol or abs(ref_start - n) < tol)):
                detected_types.append("genome_head")

            # Tail-genome: right soft-clip is genomic, alignment ends
            # near an IS copy boundary (position N or 2N)
            if (right_clip >= min_overlap
                    and aligned_len >= min_overlap
                    and (abs(ref_end - n) < tol or abs(ref_end - 2 * n) < tol)):
                detected_types.append("tail_genome")

            for junction_type in detected_types:
                if uuid not in counts:
                    counts[uuid] = {}
                counts[uuid][junction_type] = (
                    counts[uuid].get(junction_type, 0) + 1
                )

                if (junction_type == "tail_head"
                        and read.query_sequence):
                    candidate = (
                        read.mapping_quality,
                        len(read.query_sequence),
                        read.query_sequence,
                    )
                    if uuid not in th_example or candidate[:2] > th_example[uuid][:2]:
                        th_example[uuid] = candidate

                junction_reads.append({
                    "read_id": read.query_name,
                    "sample_id": sample_id,
                    "is_uuid": uuid,
                    "junction_type": junction_type,
                    "chrom": entry.chrom,
                    "start": entry.start,
                    "end": entry.end,
                    "family": entry.family,
                    "subfamily": entry.subfamily,
                    "read_length": read.query_length,
                    "alignment_start": ref_start,
                    "alignment_end": ref_end,
                    "junction_pos": junction_pos,
                    "overlap_left": junction_pos - ref_start,
                    "overlap_right": ref_end - junction_pos,
                    "mapping_quality": read.mapping_quality,
                })

        bam.close()

        # Build summary dict
        summary_by_uuid: dict[str, dict] = {}
        for uuid, entry in entries_by_uuid.items():
            uuid_counts = counts.get(uuid, {})
            summary_by_uuid[uuid] = {
                "is_uuid": uuid,
                "chrom": entry.chrom,
                "start": entry.start,
                "end": entry.end,
                "family": entry.family,
                "subfamily": entry.subfamily,
                "consensus_length": len(entry.consensus),
                "consensus": entry.consensus,
                "n_tail_head_reads": uuid_counts.get("tail_head", 0),
                "n_genome_head_reads": uuid_counts.get("genome_head", 0),
                "n_tail_genome_reads": uuid_counts.get("tail_genome", 0),
                "n_total_mapped": total_mapped.get(uuid, 0),
                "example_th_read": th_example[uuid][2] if uuid in th_example else "",
            }

        return junction_reads, summary_by_uuid

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_sample(
        self,
        tldr_table: Path,
        sample_id: str,
        fastq_path: Path,
        ref_fasta: Path | None = None,
    ) -> Optional[Path]:
        """Run circle detection for one sample.

        Args:
            tldr_table: Path to the tldr/sniffles .table.txt for this sample.
            sample_id: Sample accession.
            fastq_path: Path to the sample's FASTQ file.
            ref_fasta: Accepted for backward compatibility but ignored.

        Returns:
            Path to the sample output directory, or ``None`` on failure.
        """
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 1. Parse tldr table
        entries = self._parse_tldr_table(tldr_table)
        if not entries:
            logger.warning(
                f"No valid IS entries for {sample_id}, skipping circle detection"
            )
            return None

        entries_by_uuid = {e.uuid: e for e in entries}

        # 2. Build bait FASTA (concatemer)
        bait_path = sample_dir / f"{sample_id}_bait.fa"
        bait_path, junction_map = self._build_bait_fasta(entries, bait_path)

        # 3. Map FASTQ and collect junction reads
        bam_path = sample_dir / f"{sample_id}.circle.sorted.bam"
        bam = self._map_reads(fastq_path, bait_path, bam_path)
        if bam is None:
            logger.warning(f"Mapping failed for {sample_id}, skipping")
            return None

        junction_reads, summary_by_uuid = self._find_junction_reads(
            bam, entries_by_uuid, junction_map, sample_id,
        )

        logger.info(
            f"  {sample_id}: {len(junction_reads)} junction reads found"
        )

        # 4. Write output TSVs
        reads_tsv = sample_dir / f"{sample_id}_circle_reads.tsv"
        _write_tsv(reads_tsv, junction_reads, [
            "read_id", "sample_id", "is_uuid", "junction_type", "chrom",
            "start", "end", "family", "subfamily", "read_length",
            "alignment_start", "alignment_end", "junction_pos",
            "overlap_left", "overlap_right", "mapping_quality",
        ])

        summary_tsv = sample_dir / f"{sample_id}_circle_summary.tsv"
        summary_rows = list(summary_by_uuid.values())
        _write_tsv(summary_tsv, summary_rows, [
            "is_uuid", "chrom", "start", "end", "family", "subfamily",
            "consensus_length", "consensus", "n_tail_head_reads",
            "n_genome_head_reads", "n_tail_genome_reads", "n_total_mapped",
            "example_th_read",
        ])

        total_junctions = sum(
            s["n_tail_head_reads"] + s["n_genome_head_reads"]
            + s["n_tail_genome_reads"]
            for s in summary_by_uuid.values()
        )
        logger.info(
            f"Circle detection for {sample_id}: {total_junctions} junction "
            f"reads across {len(entries)} IS elements"
        )
        return sample_dir

    def run_batch(
        self,
        tldr_results: dict[str, Optional[Path]],
        metadata,
        fastq_dir: str | Path,
        accession_col: str = "srr_accession",
        ref_map: dict[str, Path] | None = None,
        parallel: int = 1,
    ) -> dict[str, Optional[Path]]:
        """Run circle detection for all samples.

        Args:
            tldr_results: sample_id -> tldr table path mapping.
            metadata: DataFrame with accession column.
            fastq_dir: Directory containing FASTQ files.
            accession_col: Column name for SRR accession.
            ref_map: Accepted for backward compatibility but ignored.
            parallel: Number of samples to process in parallel.

        Returns:
            Dict mapping sample_id -> circle output directory (or None).
        """
        fastq_dir = Path(fastq_dir)
        results: dict[str, Optional[Path]] = {}

        # Collect tasks
        tasks: list[tuple[str, Path, Path, Path | None]] = []
        for _, row in metadata.iterrows():
            sid = row[accession_col]
            table_path = tldr_results.get(sid)
            if not table_path:
                logger.warning(
                    f"No insertion table for {sid}, skipping circle detection"
                )
                results[sid] = None
                continue

            fq = find_fastq(fastq_dir, sid)
            if not fq:
                logger.warning(f"No FASTQ found for {sid}")
                results[sid] = None
                continue

            ref_fasta = None
            if ref_map:
                ref_fasta = ref_map.get(sid)

            tasks.append((sid, fq, table_path, ref_fasta))

        logger.info(f"Running circle detection for {len(tasks)} samples")

        if parallel <= 1:
            for sid, fq, table_path, ref_fasta in tasks:
                try:
                    out = self.run_sample(
                        tldr_table=table_path,
                        sample_id=sid,
                        fastq_path=fq,
                        ref_fasta=ref_fasta,
                    )
                    results[sid] = out
                except Exception as exc:
                    logger.error(
                        f"Circle detection failed for {sid}, skipping: {exc}"
                    )
                    results[sid] = None
        else:
            logger.info(f"Running {len(tasks)} samples with {parallel} in parallel")
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for sid, fq, table_path, ref_fasta in tasks:
                    fut = pool.submit(
                        _run_circle_worker,
                        sample_id=sid,
                        fastq_path=fq,
                        table_path=table_path,
                        ref_fasta=ref_fasta,
                        output_dir=self.output_dir,
                        threads=self.threads,
                        sort_memory=self.sort_memory,
                        min_overlap=self.min_overlap,
                        min_consensus_length=self.min_consensus_length,
                        boundary_tolerance=self.boundary_tolerance,
                        min_consensus_entropy=self.min_consensus_entropy,
                    )
                    futures[fut] = sid

                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        out = fut.result()
                        results[sid] = out
                        if out:
                            logger.info(f"  -> {out}")
                    except Exception as exc:
                        logger.error(
                            f"Circle detection worker failed for {sid}: {exc}"
                        )
                        results[sid] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(
            f"Circle detection complete: {ok}/{len(results)} samples processed"
        )
        return results


def _run_circle_worker(
    sample_id: str,
    fastq_path: Path,
    table_path: Path,
    ref_fasta: Path | None,
    output_dir: Path,
    threads: int,
    sort_memory: str,
    min_overlap: int,
    min_consensus_length: int,
    boundary_tolerance: int = DEFAULT_BOUNDARY_TOLERANCE,
    min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
    flank_length: int | None = None,  # deprecated, ignored
) -> Optional[Path]:
    """Standalone worker function for parallel circle detection."""
    finder = CircleFinder(
        output_dir=str(output_dir),
        threads=threads,
        sort_memory=sort_memory,
        min_overlap=min_overlap,
        min_consensus_length=min_consensus_length,
        boundary_tolerance=boundary_tolerance,
        min_consensus_entropy=min_consensus_entropy,
    )
    return finder.run_sample(
        tldr_table=table_path,
        sample_id=sample_id,
        fastq_path=fastq_path,
        ref_fasta=ref_fasta,
    )


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    """Write a list of dicts to a TSV file with the given column order."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"Wrote {len(rows)} rows to {path}")
