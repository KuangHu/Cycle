"""Extract IS elements with flanking regions via targeted local assembly.

After circle detection, reads mapped to IS bait sequences are extracted,
assembled per IS element with miniasm + minipolish, and the IS consensus is
located within the assembled contig to extract upstream/downstream flanks.

Output: JSON (ISExtractor-compatible) and TSV per sample.
"""

import csv
import json
import logging
import math
import re
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pysam

from .config import (
    DEFAULT_ASSEMBLY_TIMEOUT,
    DEFAULT_FLANK_SIZE,
    DEFAULT_FORMATTER_OUTPUT_DIR,
    DEFAULT_MIN_ENTROPY,
    DEFAULT_MIN_READS_FOR_ASSEMBLY,
)

logger = logging.getLogger(__name__)


class ISFormatter:
    """Extract IS elements + flanking regions from circle detection BAMs.

    For each IS element with mapped reads, performs targeted local assembly
    (miniasm + minipolish) and extracts flanking sequences from the assembled
    contig.  Falls back to the longest read when assembly fails.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_FORMATTER_OUTPUT_DIR,
        flank_size: int = DEFAULT_FLANK_SIZE,
        min_reads: int = DEFAULT_MIN_READS_FOR_ASSEMBLY,
        assembly_timeout: int = DEFAULT_ASSEMBLY_TIMEOUT,
        min_entropy: float = DEFAULT_MIN_ENTROPY,
        require_th_reads: bool = True,
        threads: int = 4,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.flank_size = flank_size
        self.require_th_reads = require_th_reads
        self.min_reads = min_reads
        self.assembly_timeout = assembly_timeout
        self.min_entropy = min_entropy
        self.threads = threads

        for tool in ("minimap2", "miniasm"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

        # minipolish is optional — we fall back to unpolished contigs
        self._have_minipolish = shutil.which("minipolish") is not None
        if not self._have_minipolish:
            logger.warning(
                "minipolish not found in PATH; will use unpolished contigs"
            )

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

    def _load_circle_summary(self, summary_tsv: Path) -> list[dict]:
        """Load circle summary TSV and return IS entries with good consensus.

        Filters out entries with no mapped reads or low-complexity consensus
        (entropy below ``min_entropy``).  When ``require_th_reads`` is True,
        also skips entries with zero tail-head junction reads.
        """
        entries: list[dict] = []
        skipped_no_reads = 0
        skipped_low_complexity = 0
        skipped_no_th = 0

        with open(summary_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                n_total = int(row.get("n_total_mapped", 0))
                if n_total == 0:
                    skipped_no_reads += 1
                    continue

                consensus = row.get("consensus", "")
                if not consensus or self._seq_entropy(consensus) < self.min_entropy:
                    skipped_low_complexity += 1
                    continue

                n_th = int(row.get("n_tail_head_reads", 0))
                if self.require_th_reads and n_th == 0:
                    skipped_no_th += 1
                    continue

                entries.append({
                    "is_uuid": row["is_uuid"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "family": row.get("family", "NA"),
                    "subfamily": row.get("subfamily", "NA"),
                    "consensus_length": int(row.get("consensus_length", 0)),
                    "consensus": consensus,
                    "n_tail_head_reads": n_th,
                    "n_genome_head_reads": int(row.get("n_genome_head_reads", 0)),
                    "n_tail_genome_reads": int(row.get("n_tail_genome_reads", 0)),
                    "n_total_mapped": n_total,
                    "example_th_read": row.get("example_th_read", ""),
                })

        total_rows = (
            len(entries) + skipped_no_reads
            + skipped_low_complexity + skipped_no_th
        )
        skip_parts = [
            f"{skipped_no_reads} no-reads",
            f"{skipped_low_complexity} low-complexity",
        ]
        if skipped_no_th:
            skip_parts.append(f"{skipped_no_th} no-TH-junction")
        logger.info(
            f"Loaded {len(entries)}/{total_rows} IS elements from "
            f"{summary_tsv.name} (skipped {', '.join(skip_parts)})"
        )
        return entries

    def _extract_all_reads(
        self, bam_paths: list[Path], uuids: set[str],
    ) -> dict[str, list[tuple[str, str, str]]]:
        """Scan all BAMs once and group reads by IS UUID.

        Each BAM is read in a single pass.  Reads are assigned to a UUID
        by parsing the bait reference name (``{uuid}__...``).  Reads are
        deduplicated by (uuid, read_name).

        Returns dict mapping uuid -> list of (read_id, sequence, qual_str).
        """
        reads_by_uuid: dict[str, list[tuple[str, str, str]]] = {
            u: [] for u in uuids
        }
        seen: dict[str, set[str]] = {u: set() for u in uuids}

        for bam_path in bam_paths:
            try:
                bam = pysam.AlignmentFile(str(bam_path), "rb")
            except Exception as exc:
                logger.warning(f"Cannot open BAM {bam_path}: {exc}")
                continue

            for read in bam.fetch():
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                ref_name = read.reference_name or ""
                # Bait headers: {uuid}__{type}__...
                sep = ref_name.find("__")
                if sep == -1:
                    continue
                uuid = ref_name[:sep]
                if uuid not in uuids:
                    continue

                seq = read.query_sequence
                if not seq:
                    continue

                name = read.query_name
                if name in seen[uuid]:
                    continue
                seen[uuid].add(name)

                qual = read.query_qualities
                if qual is not None:
                    qual_str = "".join(chr(q + 33) for q in qual)
                else:
                    qual_str = "I" * len(seq)  # default Q40

                reads_by_uuid[uuid].append((name, seq, qual_str))

            bam.close()

        total = sum(len(v) for v in reads_by_uuid.values())
        logger.info(
            f"Extracted {total} reads across {len(uuids)} IS elements "
            f"from {len(bam_paths)} BAM(s)"
        )
        return reads_by_uuid

    def _assemble_reads(
        self, fastq_path: Path, work_dir: Path,
    ) -> Optional[str]:
        """Run minimap2 overlap + miniasm + minipolish assembly pipeline.

        Returns the longest contig sequence, or None on failure.
        """
        overlaps_paf = work_dir / "overlaps.paf"
        raw_gfa = work_dir / "raw.gfa"
        polished_gfa = work_dir / "polished.gfa"

        # 1. All-vs-all overlap with minimap2
        try:
            result = subprocess.run(
                [
                    "minimap2", "-x", "ava-ont",
                    "-t", str(self.threads),
                    str(fastq_path), str(fastq_path),
                ],
                capture_output=True,
                timeout=self.assembly_timeout,
            )
            if result.returncode != 0:
                logger.debug(
                    f"minimap2 overlap failed: {result.stderr.decode().strip()}"
                )
                return None
            overlaps_paf.write_bytes(result.stdout)
        except subprocess.TimeoutExpired:
            logger.warning(f"minimap2 overlap timed out for {fastq_path.name}")
            return None
        except Exception as exc:
            logger.debug(f"minimap2 overlap error: {exc}")
            return None

        if overlaps_paf.stat().st_size == 0:
            logger.debug("No overlaps found, skipping assembly")
            return None

        # 2. miniasm — assemble unitigs
        try:
            result = subprocess.run(
                ["miniasm", "-f", str(fastq_path), str(overlaps_paf)],
                capture_output=True,
                timeout=self.assembly_timeout,
            )
            if result.returncode != 0:
                logger.debug(
                    f"miniasm failed: {result.stderr.decode().strip()}"
                )
                return None
            raw_gfa.write_bytes(result.stdout)
        except subprocess.TimeoutExpired:
            logger.warning(f"miniasm timed out for {fastq_path.name}")
            return None
        except Exception as exc:
            logger.debug(f"miniasm error: {exc}")
            return None

        contig = self._parse_gfa_longest(raw_gfa)
        if contig is None:
            return None

        # 3. minipolish — polish the assembly (optional)
        if self._have_minipolish:
            try:
                result = subprocess.run(
                    [
                        "minipolish",
                        "-t", str(self.threads),
                        str(fastq_path), str(raw_gfa),
                    ],
                    capture_output=True,
                    timeout=self.assembly_timeout * 2,
                )
                if result.returncode == 0:
                    polished_gfa.write_bytes(result.stdout)
                    polished = self._parse_gfa_longest(polished_gfa)
                    if polished is not None:
                        return polished
                    logger.debug(
                        "minipolish produced no contigs, using unpolished"
                    )
                else:
                    logger.debug(
                        f"minipolish failed: {result.stderr.decode().strip()}, "
                        f"using unpolished contig"
                    )
            except subprocess.TimeoutExpired:
                logger.warning("minipolish timed out, using unpolished contig")
            except Exception as exc:
                logger.debug(f"minipolish error: {exc}, using unpolished")

        return contig

    def _parse_gfa_longest(self, gfa_path: Path) -> Optional[str]:
        """Parse GFA S-lines and return the longest contig sequence."""
        longest = None
        longest_len = 0

        try:
            with open(gfa_path) as fh:
                for line in fh:
                    if not line.startswith("S\t"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    seq = parts[2]
                    if seq and len(seq) > longest_len:
                        longest = seq
                        longest_len = len(seq)
        except Exception as exc:
            logger.debug(f"Failed to parse GFA {gfa_path}: {exc}")
            return None

        return longest

    def _locate_is_in_contig(
        self, contig_seq: str, consensus_seq: str, work_dir: Path,
    ) -> Optional[tuple[int, int, str]]:
        """Map IS consensus to contig via minimap2 and return (start, end, strand).

        Writes temporary FASTA files for the contig and consensus, runs
        minimap2 -a, and parses the SAM output for the best alignment.
        """
        contig_fa = work_dir / "contig.fa"
        consensus_fa = work_dir / "consensus.fa"

        contig_fa.write_text(f">contig\n{contig_seq}\n")
        consensus_fa.write_text(f">consensus\n{consensus_seq}\n")

        try:
            result = subprocess.run(
                [
                    "minimap2", "-a",
                    "--secondary=no",
                    "-t", "1",
                    str(contig_fa), str(consensus_fa),
                ],
                capture_output=True,
                timeout=self.assembly_timeout,
            )
            if result.returncode != 0:
                return None
        except (subprocess.TimeoutExpired, Exception):
            return None

        # Parse SAM output
        best_start = None
        best_end = None
        best_strand = "+"
        best_mapq = -1

        for line in result.stdout.decode().splitlines():
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            if flag & 4:  # unmapped
                continue

            mapq = int(fields[4])
            pos = int(fields[3]) - 1  # SAM is 1-based
            cigar = fields[5]

            # Calculate alignment length from CIGAR
            ref_len = sum(
                int(n) for n, op in re.findall(r"(\d+)([MDNX=])", cigar)
            )

            strand = "-" if (flag & 16) else "+"

            if mapq > best_mapq:
                best_mapq = mapq
                best_start = pos
                best_end = pos + ref_len
                best_strand = strand

        if best_start is None:
            return None

        return (best_start, best_end, best_strand)

    def _extract_flanks_from_contig(
        self,
        contig_seq: str,
        is_start: int,
        is_end: int,
        flank_size: int,
    ) -> tuple[str, str]:
        """Extract upstream and downstream flanks from contig around IS position."""
        up_start = max(0, is_start - flank_size)
        upstream = contig_seq[up_start:is_start]

        down_end = min(len(contig_seq), is_end + flank_size)
        downstream = contig_seq[is_end:down_end]

        return upstream, downstream

    def _fallback_flanks_from_longest_read(
        self,
        reads: list[tuple[str, str, str]],
        consensus_seq: str,
        work_dir: Path,
        flank_size: int,
    ) -> Optional[tuple[str, str, str]]:
        """Extract IS + flanks from the longest read via alignment.

        Returns (is_sequence_from_read, upstream, downstream) or None.
        """
        if not reads:
            return None

        # Find longest read
        longest_idx = max(range(len(reads)), key=lambda i: len(reads[i][1]))
        _, longest_seq, _ = reads[longest_idx]

        # Map consensus to the longest read
        loc = self._locate_is_in_contig(longest_seq, consensus_seq, work_dir)
        if loc is None:
            return None

        is_start, is_end, strand = loc
        is_seq = longest_seq[is_start:is_end]
        upstream, downstream = self._extract_flanks_from_contig(
            longest_seq, is_start, is_end, flank_size,
        )

        return is_seq, upstream, downstream

    def _format_record(
        self,
        is_info: dict,
        is_seq: str,
        upstream: str,
        downstream: str,
        sample_id: str,
        assembly_method: str,
    ) -> dict:
        """Build ISExtractor-compatible record dict."""
        return {
            "is_id": is_info["is_uuid"],
            "sample_id": sample_id,
            "organism": "",
            "is_element": {
                "sequence": is_seq,
                "length": len(is_seq),
                "contig": is_info["chrom"],
                "start": is_info["start"],
                "end": is_info["end"],
                "strand": "+",
            },
            "flanking_upstream": {
                "sequence": upstream,
                "length": len(upstream),
            },
            "flanking_downstream": {
                "sequence": downstream,
                "length": len(downstream),
            },
            "insertion_site": {
                "contig": is_info["chrom"],
                "pos_5p": is_info["start"],
                "pos_3p": is_info["end"],
            },
            "circle_evidence": {
                "n_tail_head_reads": is_info["n_tail_head_reads"],
                "n_genome_head_reads": is_info["n_genome_head_reads"],
                "n_tail_genome_reads": is_info["n_tail_genome_reads"],
                "example_th_read": is_info.get("example_th_read", ""),
            },
            "assembly_method": assembly_method,
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_sample(
        self,
        circle_dir: Path,
        sample_id: str,
    ) -> Optional[Path]:
        """Run IS formatting for one sample's circle detection output.

        Args:
            circle_dir: Path to the sample's circle output directory
                (contains *_circle_summary.tsv and *.circle.sorted.bam).
            sample_id: Sample accession.

        Returns:
            Path to output directory containing JSON/TSV, or None on failure.
        """
        circle_dir = Path(circle_dir)

        # 1. Find and load circle summary
        summary_tsv = circle_dir / f"{sample_id}_circle_summary.tsv"
        if not summary_tsv.exists():
            # Fall back to glob for legacy layouts
            summary_files = list(circle_dir.glob("*_circle_summary.tsv"))
            if not summary_files:
                logger.warning(
                    f"No circle summary TSV found in {circle_dir}, skipping"
                )
                return None
            summary_tsv = summary_files[0]

        is_entries = self._load_circle_summary(summary_tsv)
        if not is_entries:
            logger.warning(
                f"No IS elements with reads for {sample_id}, skipping"
            )
            return None

        # 2. Discover BAMs
        bam_path = circle_dir / f"{sample_id}.circle.sorted.bam"
        if bam_path.exists():
            bam_paths = [bam_path]
        else:
            bam_paths = sorted(circle_dir.glob("*.circle.sorted.bam"))
        if not bam_paths:
            logger.warning(f"No circle BAMs found in {circle_dir}, skipping")
            return None

        logger.info(
            f"Formatting {len(is_entries)} IS elements for {sample_id} "
            f"using {len(bam_paths)} BAM(s)"
        )

        # 3. Set up output directory
        sample_out = self.output_dir / sample_id
        sample_out.mkdir(parents=True, exist_ok=True)
        assembly_dir = sample_out / "assembly"
        assembly_dir.mkdir(parents=True, exist_ok=True)

        # 4. Extract all reads from BAMs in a single pass, grouped by UUID
        all_uuids = {e["is_uuid"] for e in is_entries if e.get("consensus")}
        reads_by_uuid = self._extract_all_reads(bam_paths, all_uuids)

        # 5. Process each IS element
        records: list[dict] = []
        for is_info in is_entries:
            uuid = is_info["is_uuid"]
            uuid_prefix = uuid[:12]
            consensus_seq = is_info["consensus"]

            if not consensus_seq:
                logger.debug(f"No consensus for {uuid}, skipping")
                continue

            # Set up per-IS working directory
            is_work = assembly_dir / uuid_prefix
            is_work.mkdir(parents=True, exist_ok=True)

            # Look up pre-extracted reads for this IS
            reads = reads_by_uuid.get(uuid, [])

            is_seq = ""
            upstream = ""
            downstream = ""
            assembly_method = "none"

            if len(reads) >= self.min_reads:
                # Try assembly
                fastq_path = is_work / "reads.fq"
                with open(fastq_path, "w") as fh:
                    for read_id, seq, qual in reads:
                        fh.write(f"@{read_id}\n{seq}\n+\n{qual}\n")

                contig = self._assemble_reads(fastq_path, is_work)
                if contig is not None:
                    loc = self._locate_is_in_contig(
                        contig, consensus_seq, is_work,
                    )
                    if loc is not None:
                        is_start, is_end, strand = loc
                        is_seq = contig[is_start:is_end]
                        upstream, downstream = self._extract_flanks_from_contig(
                            contig, is_start, is_end, self.flank_size,
                        )
                        assembly_method = "assembly"
                        logger.debug(
                            f"  {uuid_prefix}: assembly OK, "
                            f"IS={len(is_seq)}bp, "
                            f"up={len(upstream)}bp, "
                            f"down={len(downstream)}bp"
                        )

            # Fallback: use longest read
            if assembly_method == "none":
                fallback = self._fallback_flanks_from_longest_read(
                    reads, consensus_seq, is_work, self.flank_size,
                )
                if fallback is not None:
                    is_seq, upstream, downstream = fallback
                    assembly_method = "longest_read"
                    logger.debug(
                        f"  {uuid_prefix}: fallback to longest read, "
                        f"IS={len(is_seq)}bp, "
                        f"up={len(upstream)}bp, "
                        f"down={len(downstream)}bp"
                    )

            record = self._format_record(
                is_info, is_seq, upstream, downstream,
                sample_id, assembly_method,
            )
            records.append(record)

        # 6. Write outputs
        if not records:
            logger.warning(f"No IS records produced for {sample_id}")
            return None

        json_path = sample_out / f"{sample_id}_is_records.json"
        with open(json_path, "w") as fh:
            json.dump(records, fh, indent=2)
        logger.info(f"Wrote {len(records)} records to {json_path}")

        tsv_path = sample_out / f"{sample_id}_is_records.tsv"
        _write_tsv(tsv_path, records)
        logger.info(f"Wrote {len(records)} rows to {tsv_path}")

        # Summary
        by_method = {}
        for r in records:
            m = r["assembly_method"]
            by_method[m] = by_method.get(m, 0) + 1
        method_str = ", ".join(f"{k}={v}" for k, v in sorted(by_method.items()))
        logger.info(
            f"IS formatting for {sample_id}: {len(records)} records ({method_str})"
        )

        return sample_out

    def run_batch(
        self,
        circle_results: dict[str, Optional[Path]],
        parallel: int = 1,
    ) -> dict[str, Optional[Path]]:
        """Run IS formatting for all samples.

        Args:
            circle_results: sample_id -> circle output directory mapping.
            parallel: Number of samples to process in parallel.

        Returns:
            Dict mapping sample_id -> formatter output directory (or None).
        """
        results: dict[str, Optional[Path]] = {}

        # Collect tasks
        tasks: list[tuple[str, Path]] = []
        for sample_id, circle_dir in circle_results.items():
            if not circle_dir:
                results[sample_id] = None
                continue
            circle_dir = Path(circle_dir)
            if not circle_dir.exists():
                logger.warning(
                    f"Circle dir does not exist for {sample_id}: {circle_dir}"
                )
                results[sample_id] = None
                continue
            tasks.append((sample_id, circle_dir))

        logger.info(f"Running IS formatting for {len(tasks)} samples")

        if parallel <= 1:
            for sample_id, circle_dir in tasks:
                try:
                    out = self.run_sample(
                        circle_dir=circle_dir,
                        sample_id=sample_id,
                    )
                    results[sample_id] = out
                except Exception as exc:
                    logger.error(
                        f"IS formatting failed for {sample_id}, skipping: {exc}"
                    )
                    results[sample_id] = None
        else:
            logger.info(
                f"Running {len(tasks)} samples with {parallel} in parallel"
            )
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for sample_id, circle_dir in tasks:
                    fut = pool.submit(
                        _run_formatter_worker,
                        sample_id=sample_id,
                        circle_dir=circle_dir,
                        output_dir=self.output_dir,
                        flank_size=self.flank_size,
                        min_reads=self.min_reads,
                        assembly_timeout=self.assembly_timeout,
                        min_entropy=self.min_entropy,
                        require_th_reads=self.require_th_reads,
                        threads=self.threads,
                    )
                    futures[fut] = sample_id

                for fut in as_completed(futures):
                    sample_id = futures[fut]
                    try:
                        out = fut.result()
                        results[sample_id] = out
                        if out:
                            logger.info(f"  -> {out}")
                    except Exception as exc:
                        logger.error(
                            f"IS formatter worker failed for {sample_id}: {exc}"
                        )
                        results[sample_id] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(
            f"IS formatting complete: {ok}/{len(results)} samples processed"
        )
        return results


# ------------------------------------------------------------------
# Module-level helpers
# ------------------------------------------------------------------

def _run_formatter_worker(
    sample_id: str,
    circle_dir: Path,
    output_dir: Path,
    flank_size: int,
    min_reads: int,
    assembly_timeout: int,
    min_entropy: float,
    require_th_reads: bool,
    threads: int,
) -> Optional[Path]:
    """Standalone worker function for parallel IS formatting."""
    formatter = ISFormatter(
        output_dir=str(output_dir),
        flank_size=flank_size,
        min_reads=min_reads,
        assembly_timeout=assembly_timeout,
        min_entropy=min_entropy,
        require_th_reads=require_th_reads,
        threads=threads,
    )
    return formatter.run_sample(
        circle_dir=circle_dir,
        sample_id=sample_id,
    )


def _write_tsv(path: Path, records: list[dict]) -> None:
    """Flatten ISExtractor-compatible records to TSV."""
    columns = [
        "is_id", "sample_id", "is_sequence", "is_length",
        "contig", "start", "end", "strand",
        "upstream_seq", "upstream_len",
        "downstream_seq", "downstream_len",
        "n_th_reads", "n_gh_reads", "n_tg_reads",
        "example_th_seq",
        "assembly_method",
    ]

    rows: list[dict] = []
    for rec in records:
        rows.append({
            "is_id": rec["is_id"],
            "sample_id": rec["sample_id"],
            "is_sequence": rec["is_element"]["sequence"],
            "is_length": rec["is_element"]["length"],
            "contig": rec["is_element"]["contig"],
            "start": rec["is_element"]["start"],
            "end": rec["is_element"]["end"],
            "strand": rec["is_element"]["strand"],
            "upstream_seq": rec["flanking_upstream"]["sequence"],
            "upstream_len": rec["flanking_upstream"]["length"],
            "downstream_seq": rec["flanking_downstream"]["sequence"],
            "downstream_len": rec["flanking_downstream"]["length"],
            "n_th_reads": rec["circle_evidence"]["n_tail_head_reads"],
            "n_gh_reads": rec["circle_evidence"]["n_genome_head_reads"],
            "n_tg_reads": rec["circle_evidence"]["n_tail_genome_reads"],
            "example_th_seq": rec["circle_evidence"]["example_th_read"],
            "assembly_method": rec["assembly_method"],
        })

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
