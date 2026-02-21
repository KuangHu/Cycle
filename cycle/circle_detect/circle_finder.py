"""Detect IS circular intermediates and insertion boundaries via bait mapping.

For each IS consensus from tldr, two bait sequences are built:

1. **Concatemer bait** ``[IS][IS]`` — reads from circular intermediates map
   across the tail-head junction (center seam).
2. **Insertion-context bait** ``[upstream_flank][IS][downstream_flank]`` — reads
   spanning the genome-head or tail-genome junctions confirm genomic integration.

Both baits are combined into a single FASTA per sample so each FASTQ is
mapped only once.
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
    DEFAULT_CIRCLE_OUTPUT_DIR,
    DEFAULT_FLANK_LENGTH,
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

    For each sample, bait sequences are built for all IS elements from
    tldr/sniffles: concatemer baits for tail-head junctions and (when a
    reference genome is provided) insertion-context baits for genome-head and
    tail-genome junctions.  The sample's FASTQ is mapped once against the
    combined bait FASTA.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_CIRCLE_OUTPUT_DIR,
        min_overlap: int = DEFAULT_MIN_JUNCTION_OVERLAP,
        min_consensus_length: int = DEFAULT_MIN_CONSENSUS_LENGTH,
        flank_length: int = DEFAULT_FLANK_LENGTH,
        min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
        threads: int = 8,
        sort_memory: str = "4G",
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_overlap = min_overlap
        self.min_consensus_length = min_consensus_length
        self.flank_length = flank_length
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

    def _extract_flanking_sequence(
        self, entry: ISEntry, fa: pysam.FastaFile,
    ) -> tuple[str, str] | None:
        """Extract upstream and downstream genomic flanks around an IS insertion.

        Uses ``pysam.FastaFile`` to fetch sequence from the reference genome.
        Handles edge cases: IS near contig boundary (clamp), contig not in
        reference (skip), flank shorter than ``min_overlap`` (skip).

        Returns:
            ``(upstream_seq, downstream_seq)`` or ``None`` if unusable.
        """
        if entry.chrom not in fa.references:
            logger.debug(
                f"Contig {entry.chrom} not in reference, "
                f"skipping IC bait for {entry.uuid}"
            )
            return None

        contig_length = fa.get_reference_length(entry.chrom)

        start = min(entry.start, entry.end)
        end = max(entry.start, entry.end)
        start = max(0, min(start, contig_length))
        end = max(0, min(end, contig_length))

        up_start = max(0, start - self.flank_length)
        down_end = min(contig_length, end + self.flank_length)

        upstream = fa.fetch(entry.chrom, up_start, start)
        downstream = fa.fetch(entry.chrom, end, down_end)

        if len(upstream) < self.min_overlap or len(downstream) < self.min_overlap:
            logger.debug(
                f"Flanks too short for {entry.uuid}: "
                f"upstream={len(upstream)}, downstream={len(downstream)}"
            )
            return None

        return upstream, downstream

    def _build_bait_fasta(
        self,
        entries: list[ISEntry],
        output_path: Path,
        ref_fasta_path: Path | None = None,
    ) -> tuple[Path, dict[str, list[dict]]]:
        """Write combined bait FASTA (concatemer + insertion-context).

        Header formats::

            >{uuid}__th__j{N}                        # tail-head bait
            >{uuid}__ic__j{J1}_j{J2}__fl{FLANK}     # insertion-context bait

        Returns:
            ``(fasta_path, junction_map)`` where *junction_map* maps
            ref_name to a list of ``{type, position, uuid}`` dicts.
        """
        output_path = Path(output_path)
        if output_path.exists():
            logger.info(f"Bait FASTA exists: {output_path}")
            return output_path, self._parse_bait_headers(output_path)

        junction_map: dict[str, list[dict]] = {}

        fa = None
        if ref_fasta_path and Path(ref_fasta_path).exists():
            try:
                fa = pysam.FastaFile(str(ref_fasta_path))
            except Exception as exc:
                logger.warning(
                    f"Could not open reference FASTA {ref_fasta_path}: {exc}"
                )

        with open(output_path, "w") as fh:
            for entry in entries:
                n = len(entry.consensus)

                # Tail-head bait (concatemer)
                th_name = f"{entry.uuid}__th__j{n}"
                fh.write(f">{th_name}\n")
                fh.write(f"{entry.consensus}{entry.consensus}\n")
                junction_map[th_name] = [
                    {"type": "tail_head", "position": n, "uuid": entry.uuid},
                ]

                # Insertion-context bait (when reference available)
                if fa is not None:
                    flanks = self._extract_flanking_sequence(entry, fa)
                    if flanks is not None:
                        upstream, downstream = flanks
                        j1 = len(upstream)
                        j2 = len(upstream) + n
                        ic_name = (
                            f"{entry.uuid}__ic__j{j1}_j{j2}"
                            f"__fl{self.flank_length}"
                        )
                        fh.write(f">{ic_name}\n")
                        fh.write(f"{upstream}{entry.consensus}{downstream}\n")
                        junction_map[ic_name] = [
                            {"type": "genome_head", "position": j1, "uuid": entry.uuid},
                            {"type": "tail_genome", "position": j2, "uuid": entry.uuid},
                        ]

        if fa is not None:
            fa.close()

        n_th = sum(1 for k in junction_map if "__th__" in k)
        n_ic = sum(1 for k in junction_map if "__ic__" in k)
        logger.info(
            f"Wrote bait FASTA with {n_th} tail-head + {n_ic} insertion-context "
            f"sequences to {output_path}"
        )
        return output_path, junction_map

    def _parse_bait_headers(self, fasta_path: Path) -> dict[str, list[dict]]:
        """Reconstruct junction_map from existing bait FASTA headers.

        Handles current ``__th__`` / ``__ic__`` formats as well as the legacy
        ``__len{N}`` format for backward compatibility.
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

                elif "__ic__j" in ref_name:
                    uuid = ref_name.split("__ic__")[0]
                    jpart = ref_name.split("__ic__j")[1].split("__fl")[0]
                    j1_str, j2_str = jpart.split("_j")
                    j1, j2 = int(j1_str), int(j2_str)
                    junction_map[ref_name] = [
                        {"type": "genome_head", "position": j1, "uuid": uuid},
                        {"type": "tail_genome", "position": j2, "uuid": uuid},
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

        Checks all junction types (tail_head, genome_head, tail_genome) for
        each aligned read.  A read is junction-spanning if its alignment starts
        at or before ``junction - min_overlap`` and ends at or after
        ``junction + min_overlap``.

        Returns:
            ``(junction_reads, summary_by_uuid)`` where *junction_reads* is a
            list of per-read dicts, and *summary_by_uuid* maps uuid to counts
            split by junction type.
        """
        junction_reads: list[dict] = []
        total_mapped: dict[str, int] = {}
        counts: dict[str, dict[str, int]] = {}  # uuid -> {type -> count}
        th_example: dict[str, str] = {}  # uuid -> first tail-head read sequence

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

            for junc in junctions:
                junction_pos = junc["position"]
                junction_type = junc["type"]

                overlap_left = junction_pos - ref_start
                overlap_right = ref_end - junction_pos

                if (ref_start <= junction_pos - self.min_overlap
                        and ref_end >= junction_pos + self.min_overlap):
                    if uuid not in counts:
                        counts[uuid] = {}
                    counts[uuid][junction_type] = (
                        counts[uuid].get(junction_type, 0) + 1
                    )

                    if (junction_type == "tail_head"
                            and uuid not in th_example
                            and read.query_sequence):
                        th_example[uuid] = read.query_sequence

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
                        "overlap_left": overlap_left,
                        "overlap_right": overlap_right,
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
                "example_th_read": th_example.get(uuid, ""),
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
            ref_fasta: Optional reference genome FASTA for insertion-context
                baits.  When ``None``, only tail-head junctions are detected.

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

        # 2. Build bait FASTA (concatemer + insertion-context)
        bait_path = sample_dir / f"{sample_id}_bait.fa"
        bait_path, junction_map = self._build_bait_fasta(
            entries, bait_path, ref_fasta,
        )

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
            ref_map: Optional sample_id -> ref_fasta path.  When provided,
                insertion-context baits are built for genome-head and
                tail-genome junction detection.
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
                        flank_length=self.flank_length,
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
    flank_length: int,
    min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
) -> Optional[Path]:
    """Standalone worker function for parallel circle detection."""
    finder = CircleFinder(
        output_dir=str(output_dir),
        threads=threads,
        sort_memory=sort_memory,
        min_overlap=min_overlap,
        min_consensus_length=min_consensus_length,
        flank_length=flank_length,
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
