"""Resolve organism names to NCBI reference genomes via the datasets CLI."""

import json
import logging
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import Optional

import pandas as pd

from .config import DEFAULT_REFERENCE_DIR

logger = logging.getLogger(__name__)


class ReferenceGenomeResolver:
    """Download one reference genome per unique organism in a metadata table.

    Resolution strategy (via ``datasets summary genome taxon``):
      1. Exact organism name with ``--reference``
      2. Genus-only (strip strain / sp. / subsp.)
      3. Log warning + skip (or use manual override TSV)

    Assembly preference: Complete Genome > Scaffold, RefSeq (GCF_) preferred.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_REFERENCE_DIR,
        override_tsv: Optional[Path] = None,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.overrides: dict[str, str] = {}
        if override_tsv and Path(override_tsv).exists():
            df = pd.read_csv(override_tsv, sep="\t")
            self.overrides = dict(zip(df["organism"], df["accession"]))
            logger.info(f"Loaded {len(self.overrides)} manual overrides")

        for tool in ("datasets", "samtools"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def resolve_all(
        self,
        metadata: pd.DataFrame,
        organism_col: str = "organism",
    ) -> dict[str, dict]:
        """Resolve and download reference genomes for all unique organisms.

        Returns:
            Dict mapping organism name → {accession, fasta, fai} or None on
            failure.
        """
        organisms = sorted(metadata[organism_col].dropna().unique())
        logger.info(
            f"{len(organisms)} unique organisms from "
            f"{len(metadata)} samples"
        )

        results: dict[str, dict] = {}
        for org in organisms:
            info = self.resolve_one(org)
            results[org] = info
            if info:
                logger.info(f"  {org} → {info['accession']}")
            else:
                logger.warning(f"  {org} → FAILED")

        ok = sum(1 for v in results.values() if v)
        logger.info(f"Resolved {ok}/{len(organisms)} organisms")
        return results

    def resolve_one(self, organism: str) -> Optional[dict]:
        """Resolve a single organism to a downloaded reference genome.

        Returns:
            Dict with keys {accession, fasta, fai} or None.
        """
        # Manual override
        if organism in self.overrides:
            accession = self.overrides[organism]
            logger.info(f"Using manual override for {organism}: {accession}")
            return self._download(accession)

        # Tier 1: exact name + --reference
        accession = self._query_taxon(organism, reference=True)
        if accession:
            return self._download(accession)

        # Tier 2: genus-only
        genus = self._extract_genus(organism)
        if genus != organism:
            logger.info(f"  Falling back to genus: {genus}")
            accession = self._query_taxon(genus, reference=True)
            if accession:
                return self._download(accession)

        logger.warning(f"No reference genome found for: {organism}")
        return None

    # ------------------------------------------------------------------
    # Internal: NCBI datasets queries
    # ------------------------------------------------------------------

    def _query_taxon(
        self,
        taxon: str,
        reference: bool = True,
    ) -> Optional[str]:
        """Query ``datasets summary genome taxon`` and pick best assembly."""
        cmd = [
            "datasets", "summary", "genome", "taxon", taxon,
            "--as-json-lines",
        ]
        if reference:
            cmd.append("--reference")

        try:
            ret = subprocess.run(
                cmd, capture_output=True, text=True, timeout=60,
            )
        except subprocess.TimeoutExpired:
            logger.warning(f"datasets timeout for: {taxon}")
            return None

        if ret.returncode != 0:
            logger.debug(f"datasets failed for {taxon}: {ret.stderr.strip()}")
            # Retry without --reference
            if reference:
                return self._query_taxon(taxon, reference=False)
            return None

        assemblies = self._parse_assemblies(ret.stdout)
        if not assemblies:
            if reference:
                return self._query_taxon(taxon, reference=False)
            return None

        best = self._pick_best(assemblies)
        return best

    @staticmethod
    def _parse_assemblies(jsonl_output: str) -> list[dict]:
        """Parse JSON-lines output from datasets summary."""
        assemblies = []
        for line in jsonl_output.strip().splitlines():
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue

            # datasets v16+ wraps assemblies in a "reports" key
            if "reports" in record:
                for r in record["reports"]:
                    assemblies.append(r)
            elif "accession" in record:
                assemblies.append(record)

        return assemblies

    @staticmethod
    def _pick_best(assemblies: list[dict]) -> Optional[str]:
        """Pick the best assembly accession from a list.

        Preference: Complete Genome > Chromosome > Scaffold > Contig,
        RefSeq (GCF_) over GenBank (GCA_).
        """
        level_rank = {
            "Complete Genome": 0,
            "Chromosome": 1,
            "Scaffold": 2,
            "Contig": 3,
        }

        def sort_key(a: dict) -> tuple:
            acc = a.get("accession", "")
            info = a.get("assembly_info", {})
            level = info.get("assembly_level", "Contig")
            is_refseq = 0 if acc.startswith("GCF_") else 1
            return (level_rank.get(level, 9), is_refseq, acc)

        assemblies.sort(key=sort_key)
        if assemblies:
            return assemblies[0].get("accession")
        return None

    @staticmethod
    def _extract_genus(organism: str) -> str:
        """Extract genus from organism name, stripping strain/sp./subsp."""
        # Remove "sp.", "subsp.", strain identifiers
        cleaned = re.sub(
            r"\s+(sp\.|subsp\.|serovar|str\.|strain)\s+.*", "", organism,
        )
        parts = cleaned.split()
        if parts:
            return parts[0]
        return organism

    # ------------------------------------------------------------------
    # Internal: download & index
    # ------------------------------------------------------------------

    def _download(self, accession: str) -> Optional[dict]:
        """Download genome via ``datasets download`` and index with samtools."""
        acc_dir = self.output_dir / accession
        fasta = acc_dir / f"{accession}_genomic.fna"
        fai = fasta.with_suffix(".fna.fai")

        # Skip if already present
        if fasta.exists() and fai.exists():
            logger.info(f"  Already downloaded: {accession}")
            return {"accession": accession, "fasta": fasta, "fai": fai}

        acc_dir.mkdir(parents=True, exist_ok=True)
        zip_path = acc_dir / f"{accession}.zip"

        # Download
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(zip_path),
        ]
        ret = subprocess.run(cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            logger.error(
                f"datasets download failed for {accession}: "
                f"{ret.stderr.strip()}"
            )
            return None

        # Unzip — find the genomic FASTA inside
        try:
            with zipfile.ZipFile(zip_path) as zf:
                fna_files = [
                    n for n in zf.namelist()
                    if n.endswith(".fna") and "genomic" in n.lower()
                ]
                if not fna_files:
                    # Fall back: any .fna file
                    fna_files = [n for n in zf.namelist() if n.endswith(".fna")]

                if not fna_files:
                    logger.error(f"No .fna in zip for {accession}")
                    return None

                # Extract the first matching .fna
                src = fna_files[0]
                with zf.open(src) as src_f, open(fasta, "wb") as dst_f:
                    dst_f.write(src_f.read())
        except zipfile.BadZipFile:
            logger.error(f"Bad zip file for {accession}")
            return None
        finally:
            zip_path.unlink(missing_ok=True)

        # Index with samtools faidx
        ret = subprocess.run(
            ["samtools", "faidx", str(fasta)],
            capture_output=True, text=True,
        )
        if ret.returncode != 0:
            logger.error(f"samtools faidx failed: {ret.stderr.strip()}")
            return None

        logger.info(f"  Downloaded and indexed {accession}")
        return {"accession": accession, "fasta": fasta, "fai": fai}
