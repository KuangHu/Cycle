"""Run Sniffles2 per sample to detect insertions as IS candidates."""

import logging
import math
import shutil
import subprocess
import uuid
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

from .config import DEFAULT_ALIGNMENT_DIR

logger = logging.getLogger(__name__)


class SnifflesRunner:
    """Run Sniffles2 once per sample to detect insertions.

    Sniffles2 is a structural variant caller that can detect insertions
    without requiring a transposon reference library. We use it as a
    faster alternative to tldr for finding IS insertion candidates.
    """

    def __init__(
        self,
        output_dir: str,
        alignment_dir: str = DEFAULT_ALIGNMENT_DIR,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.alignment_dir = Path(alignment_dir)

        if not shutil.which("sniffles"):
            raise RuntimeError("sniffles not found in PATH")

    @staticmethod
    def _seq_entropy(seq: str) -> float:
        """Shannon entropy in bits per base (0 = homopolymer, 2 = random)."""
        if not seq:
            return 0.0
        seq = seq.upper()
        n = len(seq)
        counts = Counter(seq)
        return -sum((c / n) * math.log2(c / n) for c in counts.values())

    def run_sample(
        self,
        bam: Path,
        ref_fasta: Path,
        sample_id: str,
        min_size: int = 500,
        max_size: int = 20000,
        min_support: int = 3,
        disable_qc: bool = False,
    ) -> Optional[Path]:
        """Run Sniffles2 on a single BAM for one sample.

        Args:
            bam: Sorted BAM file.
            ref_fasta: Reference genome FASTA used for alignment.
            sample_id: Sample accession (used for output directory naming).
            min_size: Minimum insertion size in bp.
            max_size: Maximum insertion size in bp.
            min_support: Minimum number of supporting reads.
            disable_qc: Disable Sniffles2 quality control filters for max sensitivity.

        Returns:
            Path to the Sniffles2-derived insertion table, or None on failure.
        """
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        output_prefix = sample_dir / sample_id
        table_path = Path(f"{output_prefix}.table.txt")
        vcf_path = Path(f"{output_prefix}.vcf")

        if table_path.exists():
            with open(table_path) as f:
                line_count = sum(1 for _ in f)
            if line_count > 1:
                logger.info(f"Sniffles2 output exists for {sample_id} ({line_count - 1} insertions): {table_path}")
                return table_path
            else:
                logger.info(f"Reprocessing {sample_id} — table exists but has no insertions")

        logger.info(f"Running Sniffles2 for {sample_id}: {bam.name}")

        # Skip if already processed with actual variants
        skip_existing = False
        if vcf_path.exists():
            with open(vcf_path) as f:
                variant_count = sum(1 for line in f if not line.startswith('#'))
            if variant_count > 0:
                logger.info(f"  {sample_id} - using existing VCF ({variant_count} variants)")
                skip_existing = True
            else:
                logger.info(f"  {sample_id} - reprocessing empty VCF")

        if not skip_existing:
            cmd = [
                "sniffles",
                "--input", str(bam),
                "--vcf", str(vcf_path),
                "--reference", str(ref_fasta),
                "--threads", "1",
                "--minsvlen", str(min_size),
                "--minsupport", str(min_support),
                "--output-rnames",
                "--allow-overwrite",
            ]
            if disable_qc:
                cmd.append("--no-qc")

            try:
                ret = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=600,
                    cwd=str(sample_dir),
                )
                if ret.returncode != 0:
                    err_msg = ret.stderr.strip() or ret.stdout.strip() or "(no output)"
                    logger.warning(f"Sniffles2 failed for {sample_id} (exit={ret.returncode}): {err_msg[:300]}")
                    return None
            except subprocess.TimeoutExpired:
                logger.warning(f"Sniffles2 timed out for {sample_id}")
                return None

        # Parse VCF and collect insertions
        insertions = []
        if vcf_path.exists():
            try:
                insertions = self._parse_vcf(vcf_path, min_size, max_size, sample_id)
            except Exception as e:
                logger.warning(f"Failed to parse VCF for {sample_id}: {e}")

        # Assemble consensus for symbolic ALTs with enough support
        assembly_dir = sample_dir / "assembly"
        symbolic = [ins for ins in insertions
                    if not ins["Sequence"] and len(ins.get("Rnames", [])) >= 3]
        n_assembled = 0
        for ins in symbolic:
            work = assembly_dir / ins["UUID"][:8]
            contig = self._assemble_insertion(
                bam_path=bam, chrom=ins["Chrom"], pos=ins["Start"],
                svlen=ins["End"] - ins["Start"],
                rnames=ins["Rnames"], work_dir=work,
            )
            if contig and self._seq_entropy(contig) >= 1.5 and len(contig) >= min_size:
                ins["Sequence"] = contig
                n_assembled += 1
        if symbolic:
            logger.info(f"  {sample_id}: assembled {n_assembled}/{len(symbolic)} symbolic ALTs")

        # Remove entries that still have no usable sequence
        insertions = [ins for ins in insertions if ins["Sequence"]]

        # Drop Rnames before writing (not needed downstream)
        for ins in insertions:
            ins.pop("Rnames", None)

        # Cleanup assembly temp files
        if assembly_dir.exists():
            shutil.rmtree(assembly_dir, ignore_errors=True)

        if insertions:
            self._write_table(insertions, table_path)
            logger.info(f"  -> {table_path} ({len(insertions)} insertions)")
            return table_path
        else:
            logger.warning(f"No insertions found for {sample_id}")
            pd.DataFrame(columns=["UUID", "Chrom", "Start", "End", "Consensus", "Support"]).to_csv(
                table_path, sep="\t", index=False
            )
            return table_path

    def _parse_vcf(
        self, vcf_path: Path, min_size: int, max_size: int, sample_id: str,
        min_entropy: float = 1.5,
    ) -> list[dict]:
        """Parse Sniffles2 VCF and extract insertions.

        Symbolic ALTs (``<INS>``) get ``Sequence=""`` and supporting read
        names stored in ``Rnames`` for downstream assembly.  Resolved
        sequences with Shannon entropy below *min_entropy* are discarded.
        """
        insertions = []
        n_symbolic = 0
        n_low_entropy = 0
        n_total = 0

        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                # Parse INFO field
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, val = item.split("=", 1)
                        info_dict[key] = val

                svtype = info_dict.get("SVTYPE", "")
                if svtype not in ("INS", "DUP", "BND"):
                    continue

                svlen = info_dict.get("SVLEN")
                if svlen:
                    svlen = abs(int(svlen))
                    if svlen < min_size or svlen > max_size:
                        continue
                elif svtype == "BND":
                    svlen = min_size
                else:
                    if alt.startswith("<"):
                        continue
                    svlen = len(alt) - len(ref)
                    if svlen < min_size or svlen > max_size:
                        continue

                n_total += 1
                support = int(info_dict.get("SUPPORT", "0"))

                # Extract read names (from --output-rnames)
                rnames_str = info_dict.get("RNAMES", "")
                rnames = [r for r in rnames_str.split(",") if r] if rnames_str else []

                if alt.startswith("<"):
                    # Symbolic ALT — mark for assembly instead of filling N's
                    n_symbolic += 1
                    sequence = ""
                else:
                    sequence = alt[len(ref):]
                    if self._seq_entropy(sequence) < min_entropy:
                        n_low_entropy += 1
                        continue

                insertions.append({
                    "UUID": str(uuid.uuid4()),
                    "Chrom": chrom,
                    "Start": pos,
                    "End": pos + svlen,
                    "Sequence": sequence,
                    "Support": support,
                    "Family": f"Sniffles2_{svtype}",
                    "Subfamily": f"{sample_id}_{chrom}_{pos}",
                    "Rnames": rnames,
                })

        logger.info(
            f"  {sample_id} VCF: {n_total} calls, {n_symbolic} symbolic, "
            f"{n_low_entropy} low-entropy filtered, {len(insertions)} kept"
        )
        return insertions

    def _assemble_insertion(
        self,
        bam_path: Path,
        chrom: str,
        pos: int,
        svlen: int,
        rnames: list[str],
        work_dir: Path,
        timeout: int = 120,
    ) -> Optional[str]:
        """Assemble a consensus for a symbolic ALT from supporting reads.

        Extracts reads by name from the BAM region around the insertion,
        then runs minimap2 overlap + miniasm to produce a consensus contig.

        Returns the longest contig sequence, or None on failure.
        """
        if not rnames:
            return None

        rname_set = set(rnames)
        work_dir.mkdir(parents=True, exist_ok=True)
        fastq_path = work_dir / "reads.fastq"

        # Fetch supporting reads from the BAM around the insertion site
        pad = max(svlen, 2000)
        start = max(0, pos - pad)
        end = pos + pad
        n_written = 0

        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                with open(fastq_path, "w") as fq:
                    for read in bam.fetch(chrom, start, end):
                        if read.query_name in rname_set and read.query_sequence:
                            qual = read.query_qualities
                            qual_str = (
                                "".join(chr(q + 33) for q in qual)
                                if qual is not None
                                else "I" * len(read.query_sequence)
                            )
                            fq.write(
                                f"@{read.query_name}\n"
                                f"{read.query_sequence}\n"
                                f"+\n"
                                f"{qual_str}\n"
                            )
                            n_written += 1
        except Exception as exc:
            logger.debug(f"Failed to extract reads at {chrom}:{pos}: {exc}")
            return None

        if n_written < 3:
            return None

        # minimap2 all-vs-all overlap
        overlaps_paf = work_dir / "overlaps.paf"
        try:
            result = subprocess.run(
                ["minimap2", "-x", "ava-ont", "-t", "1",
                 str(fastq_path), str(fastq_path)],
                capture_output=True, timeout=timeout,
            )
            if result.returncode != 0:
                return None
            overlaps_paf.write_bytes(result.stdout)
        except (subprocess.TimeoutExpired, Exception):
            return None

        if overlaps_paf.stat().st_size == 0:
            return None

        # miniasm assembly
        raw_gfa = work_dir / "raw.gfa"
        try:
            result = subprocess.run(
                ["miniasm", "-f", str(fastq_path), str(overlaps_paf)],
                capture_output=True, timeout=timeout,
            )
            if result.returncode != 0:
                return None
            raw_gfa.write_bytes(result.stdout)
        except (subprocess.TimeoutExpired, Exception):
            return None

        return self._parse_gfa_longest(raw_gfa)

    @staticmethod
    def _parse_gfa_longest(gfa_path: Path) -> Optional[str]:
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

    def _write_table(self, insertions: list[dict], output_path: Path):
        """Write insertion table in tldr-compatible format."""
        df = pd.DataFrame(insertions)
        df = df.rename(columns={"Sequence": "Consensus"})
        cols = ["UUID", "Chrom", "Start", "End", "Family", "Subfamily", "Consensus", "Support"]
        df = df[cols]
        df.to_csv(output_path, sep="\t", index=False)

    def run_batch(
        self,
        metadata: pd.DataFrame,
        ref_map: dict[str, Path],
        parallel: int = 1,
        min_size: int = 500,
        max_size: int = 20000,
        disable_qc: bool = False,
        accession_col: str = "srr_accession",
    ) -> dict[str, Optional[Path]]:
        """Run Sniffles2 for each sample.

        Args:
            metadata: DataFrame with accession column.
            ref_map: sample_id -> ref_fasta path.
            parallel: Number of samples to run in parallel.
            min_size: Minimum insertion size in bp.
            max_size: Maximum insertion size in bp.
            disable_qc: Disable Sniffles2 quality control filters.
            accession_col: Column name for SRR accession.

        Returns:
            Dict mapping sample_id -> insertion table path (or None).
        """
        results: dict[str, Optional[Path]] = {}

        # Build task list
        tasks: list[tuple[str, Path, Path]] = []
        for _, row in metadata.iterrows():
            sid = row[accession_col]
            ref_fasta = ref_map.get(sid)
            if not ref_fasta:
                logger.warning(f"No reference for {sid}, skipping Sniffles2")
                results[sid] = None
                continue

            bam = self.alignment_dir / f"{sid}.sorted.bam"
            if not bam.exists():
                logger.warning(f"BAM not found for {sid}: {bam}")
                results[sid] = None
                continue

            tasks.append((sid, bam, ref_fasta))

        logger.info(f"Running Sniffles2 for {len(tasks)} samples")

        if parallel <= 1:
            for sid, bam, ref_fasta in tasks:
                table = self.run_sample(
                    bam=bam, ref_fasta=ref_fasta,
                    sample_id=sid, min_size=min_size, max_size=max_size,
                    disable_qc=disable_qc,
                )
                results[sid] = table
        else:
            logger.info(f"Running {len(tasks)} samples with {parallel} in parallel")
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for sid, bam, ref_fasta in tasks:
                    fut = pool.submit(
                        _run_sniffles_worker,
                        bam=bam,
                        ref_fasta=ref_fasta,
                        sample_id=sid,
                        min_size=min_size,
                        max_size=max_size,
                        disable_qc=disable_qc,
                        output_dir=self.output_dir,
                        alignment_dir=self.alignment_dir,
                    )
                    futures[fut] = sid

                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        table = fut.result()
                        results[sid] = table
                        if table:
                            logger.info(f"  -> {table}")
                    except Exception as e:
                        logger.error(f"Sniffles2 worker failed for {sid}: {e}")
                        results[sid] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(f"Sniffles2 complete: {ok}/{len(results)} samples succeeded")
        return results


def _run_sniffles_worker(
    bam: Path,
    ref_fasta: Path,
    sample_id: str,
    min_size: int,
    max_size: int,
    disable_qc: bool,
    output_dir: Path,
    alignment_dir: Path,
) -> Optional[Path]:
    """Standalone worker function for parallel Sniffles2 execution."""
    runner = SnifflesRunner(
        output_dir=str(output_dir),
        alignment_dir=str(alignment_dir),
    )
    return runner.run_sample(
        bam=bam,
        ref_fasta=ref_fasta,
        sample_id=sample_id,
        min_size=min_size,
        max_size=max_size,
        disable_qc=disable_qc,
    )
