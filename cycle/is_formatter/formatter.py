"""Extract IS elements with flanking regions via targeted local assembly.

After circle detection, reads mapped to IS bait sequences are extracted,
assembled per IS element with miniasm + minipolish, and the IS consensus is
located within the assembled contig to extract upstream/downstream flanks.

Output: JSON (ISExtractor-compatible) and TSV per organism.
"""

import csv
import json
import logging
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pysam

from ..utils import slugify
from .config import (
    DEFAULT_ASSEMBLY_TIMEOUT,
    DEFAULT_FLANK_SIZE,
    DEFAULT_FORMATTER_OUTPUT_DIR,
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
        threads: int = 4,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.flank_size = flank_size
        self.min_reads = min_reads
        self.assembly_timeout = assembly_timeout
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

    def _load_circle_summary(self, summary_tsv: Path) -> list[dict]:
        """Load circle summary TSV and return IS entries with any reads."""
        entries: list[dict] = []
        with open(summary_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                n_total = int(row.get("n_total_mapped", 0))
                if n_total == 0:
                    continue
                entries.append({
                    "is_uuid": row["is_uuid"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "family": row.get("family", "NA"),
                    "subfamily": row.get("subfamily", "NA"),
                    "consensus_length": int(row.get("consensus_length", 0)),
                    "consensus": row.get("consensus", ""),
                    "n_tail_head_reads": int(row.get("n_tail_head_reads", 0)),
                    "n_genome_head_reads": int(row.get("n_genome_head_reads", 0)),
                    "n_tail_genome_reads": int(row.get("n_tail_genome_reads", 0)),
                    "n_total_mapped": n_total,
                })
        logger.info(
            f"Loaded {len(entries)} IS elements with reads from "
            f"{summary_tsv.name}"
        )
        return entries

    def _extract_reads_for_is(
        self, bam_paths: list[Path], is_uuid: str,
    ) -> list[tuple[str, str, str]]:
        """Extract reads mapped to any bait containing *is_uuid*.

        Returns list of (read_id, sequence, quality_string) tuples,
        deduplicated by read name.
        """
        seen: set[str] = set()
        reads: list[tuple[str, str, str]] = []

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
                # Bait headers contain the uuid before the __ separator
                if not ref_name.startswith(is_uuid + "__"):
                    continue

                seq = read.query_sequence
                qual = read.query_qualities
                if not seq:
                    continue

                name = read.query_name
                if name in seen:
                    continue
                seen.add(name)

                # Convert quality array to FASTQ quality string
                if qual is not None:
                    qual_str = "".join(chr(q + 33) for q in qual)
                else:
                    qual_str = "I" * len(seq)  # default Q40

                reads.append((name, seq, qual_str))

            bam.close()

        return reads

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
        organism: str,
        assembly_method: str,
    ) -> dict:
        """Build ISExtractor-compatible record dict."""
        return {
            "is_id": is_info["is_uuid"],
            "organism": organism,
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
            },
            "assembly_method": assembly_method,
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_organism(
        self,
        circle_dir: Path,
        organism: str,
    ) -> Optional[Path]:
        """Run IS formatting for one organism's circle detection output.

        Args:
            circle_dir: Path to the organism's circle output directory
                (contains *_circle_summary.tsv and *.circle.sorted.bam).
            organism: Organism name (for output naming and record metadata).

        Returns:
            Path to output directory containing JSON/TSV, or None on failure.
        """
        circle_dir = Path(circle_dir)
        slug = slugify(organism)

        # 1. Find and load circle summary
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
                f"No IS elements with reads for {organism}, skipping"
            )
            return None

        # 2. Discover BAMs
        bam_paths = sorted(circle_dir.glob("*.circle.sorted.bam"))
        if not bam_paths:
            logger.warning(f"No circle BAMs found in {circle_dir}, skipping")
            return None

        logger.info(
            f"Formatting {len(is_entries)} IS elements for {organism} "
            f"using {len(bam_paths)} BAM(s)"
        )

        # 3. Set up output directory
        org_out = self.output_dir / slug
        org_out.mkdir(parents=True, exist_ok=True)
        assembly_dir = org_out / "assembly"
        assembly_dir.mkdir(parents=True, exist_ok=True)

        # 4. Process each IS element
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

            # Extract reads for this IS
            reads = self._extract_reads_for_is(bam_paths, uuid)
            logger.debug(f"  {uuid_prefix}: {len(reads)} reads extracted")

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
                organism, assembly_method,
            )
            records.append(record)

        # 5. Write outputs
        if not records:
            logger.warning(f"No IS records produced for {organism}")
            return None

        json_path = org_out / f"{slug}_is_records.json"
        with open(json_path, "w") as fh:
            json.dump(records, fh, indent=2)
        logger.info(f"Wrote {len(records)} records to {json_path}")

        tsv_path = org_out / f"{slug}_is_records.tsv"
        _write_tsv(tsv_path, records)
        logger.info(f"Wrote {len(records)} rows to {tsv_path}")

        # Summary
        by_method = {}
        for r in records:
            m = r["assembly_method"]
            by_method[m] = by_method.get(m, 0) + 1
        method_str = ", ".join(f"{k}={v}" for k, v in sorted(by_method.items()))
        logger.info(
            f"IS formatting for {organism}: {len(records)} records ({method_str})"
        )

        return org_out

    def run_batch(
        self,
        circle_results: dict[str, Optional[Path]],
        metadata,
        parallel: int = 1,
    ) -> dict[str, Optional[Path]]:
        """Run IS formatting for all organism groups.

        Args:
            circle_results: Organism -> circle output directory mapping.
            metadata: DataFrame with organism column (used for organism names).
            parallel: Number of organisms to process in parallel.

        Returns:
            Dict mapping organism -> formatter output directory (or None).
        """
        results: dict[str, Optional[Path]] = {}

        # Collect tasks
        tasks: list[tuple[str, Path]] = []
        for organism, circle_dir in circle_results.items():
            if not circle_dir:
                results[organism] = None
                continue
            circle_dir = Path(circle_dir)
            if not circle_dir.exists():
                logger.warning(
                    f"Circle dir does not exist for {organism}: {circle_dir}"
                )
                results[organism] = None
                continue
            tasks.append((organism, circle_dir))

        logger.info(f"Running IS formatting for {len(tasks)} organism groups")

        if parallel <= 1:
            for organism, circle_dir in tasks:
                try:
                    out = self.run_organism(
                        circle_dir=circle_dir,
                        organism=organism,
                    )
                    results[organism] = out
                except Exception as exc:
                    logger.error(
                        f"IS formatting failed for {organism}, skipping: {exc}"
                    )
                    results[organism] = None
        else:
            logger.info(
                f"Running {len(tasks)} organisms with {parallel} in parallel"
            )
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for organism, circle_dir in tasks:
                    fut = pool.submit(
                        _run_formatter_worker,
                        organism=organism,
                        circle_dir=circle_dir,
                        output_dir=self.output_dir,
                        flank_size=self.flank_size,
                        min_reads=self.min_reads,
                        assembly_timeout=self.assembly_timeout,
                        threads=self.threads,
                    )
                    futures[fut] = organism

                for fut in as_completed(futures):
                    organism = futures[fut]
                    try:
                        out = fut.result()
                        results[organism] = out
                        if out:
                            logger.info(f"  -> {out}")
                    except Exception as exc:
                        logger.error(
                            f"IS formatter worker failed for {organism}: {exc}"
                        )
                        results[organism] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(
            f"IS formatting complete: {ok}/{len(results)} organism groups "
            f"processed"
        )
        return results


# ------------------------------------------------------------------
# Module-level helpers
# ------------------------------------------------------------------

def _run_formatter_worker(
    organism: str,
    circle_dir: Path,
    output_dir: Path,
    flank_size: int,
    min_reads: int,
    assembly_timeout: int,
    threads: int,
) -> Optional[Path]:
    """Standalone worker function for parallel IS formatting."""
    formatter = ISFormatter(
        output_dir=str(output_dir),
        flank_size=flank_size,
        min_reads=min_reads,
        assembly_timeout=assembly_timeout,
        threads=threads,
    )
    return formatter.run_organism(
        circle_dir=circle_dir,
        organism=organism,
    )


def _write_tsv(path: Path, records: list[dict]) -> None:
    """Flatten ISExtractor-compatible records to TSV."""
    columns = [
        "is_id", "organism", "is_sequence", "is_length",
        "contig", "start", "end", "strand",
        "upstream_seq", "upstream_len",
        "downstream_seq", "downstream_len",
        "n_th_reads", "n_gh_reads", "n_tg_reads",
        "assembly_method",
    ]

    rows: list[dict] = []
    for rec in records:
        rows.append({
            "is_id": rec["is_id"],
            "organism": rec["organism"],
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
            "assembly_method": rec["assembly_method"],
        })

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
