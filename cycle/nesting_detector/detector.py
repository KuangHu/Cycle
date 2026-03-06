"""Detect IS-within-IS nesting events by pairwise alignment of IS110 elements.

Aligns extended sequences (upstream + IS + downstream) all-vs-all with minimap2,
then identifies pairs where a longer "host" element contains insertion(s) relative
to a shorter "core" element — evidence of IS-within-IS nesting.
"""

import csv
import json
import logging
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .config import (
    DEFAULT_FLANKING_PAD,
    DEFAULT_MIN_BLOCK_LENGTH,
    DEFAULT_MIN_IDENTITY,
    DEFAULT_MIN_INSERTION_SIZE,
    DEFAULT_MIN_LENGTH_RATIO,
    DEFAULT_MM2_PRESET,
    DEFAULT_NESTING_OUTPUT_DIR,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class AlignmentBlock:
    """One contiguous alignment block from a PAF record."""

    query_start: int
    query_end: int
    target_start: int
    target_end: int
    strand: str
    residue_matches: int
    block_length: int

    @property
    def identity(self) -> float:
        return self.residue_matches / self.block_length if self.block_length else 0.0


@dataclass
class Insertion:
    """A detected insertion in the host relative to the core."""

    host_start: int
    host_end: int
    core_position: int
    insertion_size: int
    inserted_sequence: str = ""


@dataclass
class NestingEvent:
    """A nesting event: one IS element inserted into another."""

    host_is_id: str
    core_is_id: str
    host_sample_id: str
    core_sample_id: str
    host_is_length: int
    core_is_length: int
    n_aligned_blocks: int
    total_aligned_bp: int
    alignment_coverage: float      # fraction of core covered
    mean_identity: float
    insertions: list[Insertion] = field(default_factory=list)
    aligned_blocks: list[dict] = field(default_factory=list)
    # Each block: {core_start, core_end, host_start, host_end, identity}

    @property
    def n_insertions(self) -> int:
        return len(self.insertions)

    @property
    def total_insertion_bp(self) -> int:
        return sum(ins.insertion_size for ins in self.insertions)


# ---------------------------------------------------------------------------
# Main detector class
# ---------------------------------------------------------------------------

class NestingDetector:
    """Detect IS-within-IS nesting by pairwise alignment of IS110 elements.

    Workflow
    --------
    1. Load IS110 records and build extended sequences (flanking + IS).
    2. Write FASTA and run minimap2 all-vs-all (``-X`` flag).
    3. Parse PAF output and detect split-alignment patterns indicating nesting.
    4. Extract inserted sequences and write JSON + TSV output.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_NESTING_OUTPUT_DIR,
        mm2_preset: str = DEFAULT_MM2_PRESET,
        min_identity: float = DEFAULT_MIN_IDENTITY,
        min_block_length: int = DEFAULT_MIN_BLOCK_LENGTH,
        min_insertion_size: int = DEFAULT_MIN_INSERTION_SIZE,
        flanking_pad: int = DEFAULT_FLANKING_PAD,
        min_length_ratio: float = DEFAULT_MIN_LENGTH_RATIO,
        threads: int = 4,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.mm2_preset = mm2_preset
        self.min_identity = min_identity
        self.min_block_length = min_block_length
        self.min_insertion_size = min_insertion_size
        self.flanking_pad = flanking_pad
        self.min_length_ratio = min_length_ratio
        self.threads = threads

        # Validate minimap2 is available
        if not shutil.which("minimap2"):
            raise RuntimeError("minimap2 not found in PATH")

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def run(self, records_path: Path, skip_existing: bool = True) -> Path:
        """Run the full nesting detection pipeline.

        Parameters
        ----------
        records_path : Path
            Path to ``is110_circular_records.json``.
        skip_existing : bool
            Skip if output files already exist.

        Returns
        -------
        Path to the output directory.
        """
        records_path = Path(records_path)
        events_json = self.output_dir / "nesting_events.json"

        if skip_existing and events_json.exists():
            logger.info("Output already exists, skipping: %s", events_json)
            return self.output_dir

        # 1. Load records and build extended sequences
        records = self._load_records(records_path)
        if not records:
            logger.warning("No records loaded — nothing to do")
            return self.output_dir

        sequences = self._build_extended_sequences(records)
        if len(sequences) < 2:
            logger.warning("Need at least 2 sequences for pairwise alignment")
            return self.output_dir

        # 2. Write FASTA and run minimap2
        fasta_path = self._write_fasta(sequences)
        paf_path = self._run_minimap2(fasta_path)
        if paf_path is None:
            logger.error("minimap2 failed — aborting")
            return self.output_dir

        # 3. Parse PAF and detect nesting
        blocks_by_pair = self._parse_paf(paf_path)
        events = self._detect_nesting(blocks_by_pair, sequences)

        # 4. Write output
        self._write_events_json(events, sequences)
        self._write_events_tsv(events)

        logger.info(
            "Nesting detection complete: %d events → %s",
            len(events), self.output_dir,
        )
        return self.output_dir

    # ------------------------------------------------------------------
    # Internal methods
    # ------------------------------------------------------------------

    def _load_records(self, records_path: Path) -> list[dict]:
        """Load IS110 records from JSON."""
        try:
            with open(records_path) as fh:
                records = json.load(fh)
        except Exception as exc:
            logger.error("Failed to load records from %s: %s", records_path, exc)
            return []

        logger.info("Loaded %d IS110 records from %s", len(records), records_path.name)
        return records

    def _build_extended_sequences(
        self, records: list[dict],
    ) -> dict[str, dict]:
        """Build extended sequences: upstream + IS + downstream.

        Returns dict keyed by is_id with fields:
            extended_seq, is_seq, is_length, sample_id, ext_length
        """
        sequences = {}
        skipped = 0

        for rec in records:
            is_id = rec.get("is_id", "")
            sample_id = rec.get("sample_id", "")
            is_seq = (rec.get("is_element") or {}).get("sequence", "")

            if not is_seq:
                skipped += 1
                continue

            up_seq = (rec.get("flanking_upstream") or {}).get("sequence", "")
            dn_seq = (rec.get("flanking_downstream") or {}).get("sequence", "")

            # Trim flanking to pad size
            up_pad = up_seq[-self.flanking_pad:] if up_seq else ""
            dn_pad = dn_seq[:self.flanking_pad] if dn_seq else ""

            extended = up_pad + is_seq + dn_pad

            sequences[is_id] = {
                "extended_seq": extended,
                "is_seq": is_seq,
                "is_length": len(is_seq),
                "ext_length": len(extended),
                "sample_id": sample_id,
                "up_pad_len": len(up_pad),
                "dn_pad_len": len(dn_pad),
            }

        if skipped:
            logger.warning("Skipped %d records with no IS sequence", skipped)
        logger.info("Built %d extended sequences", len(sequences))
        return sequences

    def _write_fasta(self, sequences: dict[str, dict]) -> Path:
        """Write extended sequences to FASTA."""
        fasta_path = self.output_dir / "extended_sequences.fasta"
        with open(fasta_path, "w") as fh:
            for is_id, info in sequences.items():
                fh.write(f">{is_id}\n{info['extended_seq']}\n")
        logger.info("Wrote %d sequences to %s", len(sequences), fasta_path.name)
        return fasta_path

    def _run_minimap2(self, fasta_path: Path) -> Optional[Path]:
        """Run minimap2 all-vs-all and write PAF output."""
        paf_path = self.output_dir / "all_vs_all.paf"

        cmd = [
            "minimap2",
            "-x", self.mm2_preset,
            "-X",                       # all-vs-all (skip self hits)
            "-c",                       # output CIGAR
            "--eqx",                    # use =/X in CIGAR
            "-t", str(self.threads),
            str(fasta_path),
            str(fasta_path),
        ]

        logger.info("Running minimap2 all-vs-all (%d threads)", self.threads)
        try:
            with open(paf_path, "w") as out_fh:
                proc = subprocess.run(
                    cmd,
                    stdout=out_fh,
                    stderr=subprocess.PIPE,
                    check=True,
                )
        except subprocess.CalledProcessError as exc:
            logger.error(
                "minimap2 failed (rc=%d): %s",
                exc.returncode, exc.stderr.decode().strip(),
            )
            return None

        # Count alignments
        n_lines = sum(1 for _ in open(paf_path))
        logger.info("minimap2 produced %d alignments → %s", n_lines, paf_path.name)
        return paf_path

    def _parse_paf(
        self, paf_path: Path,
    ) -> dict[tuple[str, str], list[AlignmentBlock]]:
        """Parse PAF into alignment blocks grouped by (query, target) pair.

        When a CIGAR string is present (``cg:Z:`` tag), large alignments are
        decomposed into sub-blocks by splitting at large D/I operations
        (>= *min_insertion_size*).  This lets the downstream gap-detection
        logic naturally find insertions that minimap2 represents as internal
        indels rather than separate alignment records.
        """
        blocks_by_pair: dict[tuple[str, str], list[AlignmentBlock]] = {}

        with open(paf_path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue

                query_name = fields[0]
                query_start = int(fields[2])
                strand = fields[4]
                target_name = fields[5]
                target_start = int(fields[7])
                residue_matches = int(fields[9])
                block_length = int(fields[10])

                # Skip self-alignments (shouldn't appear with -X, but safe)
                if query_name == target_name:
                    continue

                # Look for CIGAR tag
                cigar = None
                for f in fields[12:]:
                    if f.startswith("cg:Z:"):
                        cigar = f[5:]
                        break

                if cigar:
                    sub_blocks = self._decompose_cigar(
                        cigar, query_start, target_start, strand,
                    )
                else:
                    # No CIGAR — use the whole alignment as one block
                    sub_blocks = [AlignmentBlock(
                        query_start=query_start,
                        query_end=int(fields[3]),
                        target_start=target_start,
                        target_end=int(fields[8]),
                        strand=strand,
                        residue_matches=residue_matches,
                        block_length=block_length,
                    )]

                key = (query_name, target_name)
                if key not in blocks_by_pair:
                    blocks_by_pair[key] = []

                for blk in sub_blocks:
                    if blk.block_length < self.min_block_length:
                        continue
                    if blk.identity < self.min_identity:
                        continue
                    blocks_by_pair[key].append(blk)

        logger.info(
            "Parsed %d alignment pairs from PAF",
            len(blocks_by_pair),
        )
        return blocks_by_pair

    def _decompose_cigar(
        self,
        cigar: str,
        query_start: int,
        target_start: int,
        strand: str,
    ) -> list[AlignmentBlock]:
        """Split a CIGAR alignment into sub-blocks at large indels.

        Large D operations (>= *min_insertion_size*) indicate sequence present
        in the target but absent from the query; large I operations indicate
        the reverse.  Splitting at these creates sub-blocks whose inter-block
        gaps are then detected by ``_find_insertions``.
        """
        ops = re.findall(r"(\d+)([=XIDMSH])", cigar)

        q_pos = query_start
        t_pos = target_start

        # Accumulate current sub-block stats
        blk_q_start = q_pos
        blk_t_start = t_pos
        blk_matches = 0
        blk_aligned = 0  # aligned columns (=, X, M only)

        blocks: list[AlignmentBlock] = []

        def _flush_block() -> None:
            if blk_aligned > 0:
                blocks.append(AlignmentBlock(
                    query_start=blk_q_start,
                    query_end=q_pos,
                    target_start=blk_t_start,
                    target_end=t_pos,
                    strand=strand,
                    residue_matches=blk_matches,
                    block_length=blk_aligned,
                ))

        for length_str, op in ops:
            length = int(length_str)

            if op in ("=", "M"):
                blk_matches += length
                blk_aligned += length
                q_pos += length
                t_pos += length

            elif op == "X":
                blk_aligned += length
                q_pos += length
                t_pos += length

            elif op == "D":
                # Deletion in query = extra bases in target
                if length >= self.min_insertion_size:
                    _flush_block()
                    t_pos += length
                    blk_q_start = q_pos
                    blk_t_start = t_pos
                    blk_matches = 0
                    blk_aligned = 0
                else:
                    # Small deletion — keep within current sub-block
                    t_pos += length

            elif op == "I":
                # Insertion in query = extra bases in query
                if length >= self.min_insertion_size:
                    _flush_block()
                    q_pos += length
                    blk_q_start = q_pos
                    blk_t_start = t_pos
                    blk_matches = 0
                    blk_aligned = 0
                else:
                    q_pos += length

            # S/H: soft/hard clip — only at ends, skip
            elif op == "S":
                q_pos += length
            # H doesn't consume either

        # Final sub-block
        _flush_block()

        return blocks

    def _detect_nesting(
        self,
        blocks_by_pair: dict[tuple[str, str], list[AlignmentBlock]],
        sequences: dict[str, dict],
    ) -> list[NestingEvent]:
        """Detect nesting events from alignment blocks.

        For each pair, the shorter element is the "core" and the longer is
        the "host".  We look for gaps between aligned blocks on the host that
        are larger than the corresponding gaps on the core.
        """
        events = []
        seen_pairs: set[tuple[str, str]] = set()

        for (qname, tname), blocks in blocks_by_pair.items():
            # Deduplicate: only process each unordered pair once
            canonical = tuple(sorted([qname, tname]))
            if canonical in seen_pairs:
                continue
            seen_pairs.add(canonical)

            q_info = sequences.get(qname)
            t_info = sequences.get(tname)
            if not q_info or not t_info:
                continue

            # Determine host (longer) and core (shorter)
            q_ext_len = q_info["ext_length"]
            t_ext_len = t_info["ext_length"]

            if q_ext_len >= t_ext_len:
                host_id, core_id = qname, tname
                host_info, core_info = q_info, t_info
            else:
                host_id, core_id = tname, qname
                host_info, core_info = t_info, q_info

            # Check minimum length ratio
            ratio = host_info["ext_length"] / core_info["ext_length"]
            if ratio < self.min_length_ratio:
                continue

            # Collect all blocks for this pair (may come from either direction)
            all_blocks = list(blocks)
            reverse_key = (tname, qname)
            if reverse_key in blocks_by_pair:
                all_blocks.extend(blocks_by_pair[reverse_key])

            # Orient blocks so query=core, target=host
            oriented = self._orient_blocks(all_blocks, core_id, host_id,
                                           qname, tname, sequences)
            if not oriented:
                continue

            # Sort by position on core (query)
            oriented.sort(key=lambda b: b.query_start)

            # Check collinearity on host
            if not self._is_collinear(oriented):
                continue

            # Detect insertions from gaps between blocks
            insertions = self._find_insertions(oriented, host_info)

            if not insertions:
                continue

            # Compute summary stats
            total_aligned = sum(b.query_end - b.query_start for b in oriented)
            coverage = total_aligned / core_info["ext_length"]
            total_matches = sum(b.residue_matches for b in oriented)
            total_block_len = sum(b.block_length for b in oriented)
            mean_id = total_matches / total_block_len if total_block_len else 0.0

            # Store aligned blocks for visualization
            block_dicts = [
                {
                    "core_start": b.query_start,
                    "core_end": b.query_end,
                    "host_start": b.target_start,
                    "host_end": b.target_end,
                    "identity": round(b.identity, 4),
                }
                for b in oriented
            ]

            event = NestingEvent(
                host_is_id=host_id,
                core_is_id=core_id,
                host_sample_id=host_info["sample_id"],
                core_sample_id=core_info["sample_id"],
                host_is_length=host_info["is_length"],
                core_is_length=core_info["is_length"],
                n_aligned_blocks=len(oriented),
                total_aligned_bp=total_aligned,
                alignment_coverage=round(coverage, 4),
                mean_identity=round(mean_id, 4),
                insertions=insertions,
                aligned_blocks=block_dicts,
            )
            events.append(event)

        logger.info("Detected %d nesting events", len(events))
        return events

    def _orient_blocks(
        self,
        blocks: list[AlignmentBlock],
        core_id: str,
        host_id: str,
        qname: str,
        tname: str,
        sequences: dict[str, dict],
    ) -> list[AlignmentBlock]:
        """Re-orient blocks so query coords = core, target coords = host."""
        oriented = []
        for b in blocks:
            # In PAF, query is fields[0] and target is fields[5]
            # If qname == core_id, coords are already right
            if qname == core_id:
                oriented.append(AlignmentBlock(
                    query_start=b.query_start,
                    query_end=b.query_end,
                    target_start=b.target_start,
                    target_end=b.target_end,
                    strand=b.strand,
                    residue_matches=b.residue_matches,
                    block_length=b.block_length,
                ))
            else:
                # Swap query and target
                oriented.append(AlignmentBlock(
                    query_start=b.target_start,
                    query_end=b.target_end,
                    target_start=b.query_start,
                    target_end=b.query_end,
                    strand=b.strand,
                    residue_matches=b.residue_matches,
                    block_length=b.block_length,
                ))
        return oriented

    @staticmethod
    def _is_collinear(blocks: list[AlignmentBlock]) -> bool:
        """Check that blocks are in the same order on both sequences."""
        if len(blocks) < 2:
            return True
        for i in range(1, len(blocks)):
            if blocks[i].target_start < blocks[i - 1].target_start:
                return False
        return True

    def _find_insertions(
        self,
        blocks: list[AlignmentBlock],
        host_info: dict,
    ) -> list[Insertion]:
        """Find insertions from gaps between consecutive alignment blocks."""
        insertions = []
        host_seq = host_info["extended_seq"]

        for i in range(1, len(blocks)):
            prev = blocks[i - 1]
            curr = blocks[i]

            gap_on_core = curr.query_start - prev.query_end
            gap_on_host = curr.target_start - prev.target_end

            insertion_size = gap_on_host - max(gap_on_core, 0)

            if insertion_size >= self.min_insertion_size:
                host_start = prev.target_end
                host_end = curr.target_start
                inserted_seq = host_seq[host_start:host_end]

                insertions.append(Insertion(
                    host_start=host_start,
                    host_end=host_end,
                    core_position=prev.query_end,
                    insertion_size=insertion_size,
                    inserted_sequence=inserted_seq,
                ))

        return insertions

    def _write_events_json(
        self, events: list[NestingEvent], sequences: dict[str, dict],
    ) -> None:
        """Write detailed nesting events to JSON."""
        path = self.output_dir / "nesting_events.json"

        records = []
        for ev in events:
            records.append({
                "host_is_id": ev.host_is_id,
                "core_is_id": ev.core_is_id,
                "host_sample_id": ev.host_sample_id,
                "core_sample_id": ev.core_sample_id,
                "host_is_length": ev.host_is_length,
                "core_is_length": ev.core_is_length,
                "n_aligned_blocks": ev.n_aligned_blocks,
                "total_aligned_bp": ev.total_aligned_bp,
                "alignment_coverage": ev.alignment_coverage,
                "mean_identity": ev.mean_identity,
                "n_insertions": ev.n_insertions,
                "total_insertion_bp": ev.total_insertion_bp,
                "insertions": [
                    {
                        "host_start": ins.host_start,
                        "host_end": ins.host_end,
                        "core_position": ins.core_position,
                        "insertion_size": ins.insertion_size,
                        "inserted_sequence": ins.inserted_sequence,
                    }
                    for ins in ev.insertions
                ],
                "aligned_blocks": ev.aligned_blocks,
            })

        with open(path, "w") as fh:
            json.dump(records, fh, indent=2)
        logger.info("Wrote %d nesting events to %s", len(records), path.name)

    def _write_events_tsv(self, events: list[NestingEvent]) -> None:
        """Write summary TSV with one row per insertion."""
        path = self.output_dir / "nesting_events_summary.tsv"
        columns = [
            "host_is_id", "core_is_id", "host_sample", "core_sample",
            "host_len", "core_len", "insertion_pos_on_core",
            "insertion_size", "n_blocks", "coverage", "identity",
        ]

        rows = []
        for ev in events:
            for ins in ev.insertions:
                rows.append({
                    "host_is_id": ev.host_is_id,
                    "core_is_id": ev.core_is_id,
                    "host_sample": ev.host_sample_id,
                    "core_sample": ev.core_sample_id,
                    "host_len": ev.host_is_length,
                    "core_len": ev.core_is_length,
                    "insertion_pos_on_core": ins.core_position,
                    "insertion_size": ins.insertion_size,
                    "n_blocks": ev.n_aligned_blocks,
                    "coverage": ev.alignment_coverage,
                    "identity": ev.mean_identity,
                })

        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

        logger.info("Wrote %d rows to %s", len(rows), path.name)
