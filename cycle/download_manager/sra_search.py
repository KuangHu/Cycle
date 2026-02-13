"""Search NCBI SRA via Entrez E-utilities and collect run metadata."""

import time
import logging
import xml.etree.ElementTree as ET
from io import BytesIO
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import Entrez

from .config import (
    DEFAULT_SRA_QUERY,
    SRA_DB,
    BATCH_SIZE,
    REQUEST_DELAY,
)

logger = logging.getLogger(__name__)


def _text(elem, tag, default=""):
    """Get text of a child element, or default."""
    child = elem.find(tag)
    if child is not None and child.text:
        return child.text.strip()
    return default


def _attr(elem, key, default=""):
    """Get an attribute value from an element."""
    return elem.get(key, default)


class SRASearcher:
    """Search SRA and collect run-level metadata into a DataFrame."""

    def __init__(self, email: str, api_key: Optional[str] = None):
        """
        Args:
            email: Required by NCBI Entrez.
            api_key: NCBI API key (allows 10 req/s instead of 3).
        """
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        self.delay = 0.11 if api_key else REQUEST_DELAY

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def search(
        self,
        query: str = DEFAULT_SRA_QUERY,
        max_results: int = 0,
    ) -> list[str]:
        """Return list of SRA UIDs matching *query*.

        Args:
            query: Entrez search string.
            max_results: Cap on total UIDs (0 = all).
        """
        # First pass: get total count
        handle = Entrez.esearch(db=SRA_DB, term=query, retmax=0, usehistory="y")
        record = Entrez.read(handle)
        handle.close()

        total = int(record["Count"])
        web_env = record["WebEnv"]
        query_key = record["QueryKey"]
        logger.info(f"SRA search returned {total} results")

        if max_results > 0:
            total = min(total, max_results)

        # Paginate through history
        uids: list[str] = []
        for start in range(0, total, BATCH_SIZE):
            batch = min(BATCH_SIZE, total - start)
            handle = Entrez.esearch(
                db=SRA_DB,
                term=query,
                retstart=start,
                retmax=batch,
                webenv=web_env,
                query_key=query_key,
            )
            rec = Entrez.read(handle)
            handle.close()
            uids.extend(rec["IdList"])
            logger.info(f"  fetched UIDs {start}–{start + len(rec['IdList'])}")
            time.sleep(self.delay)

        logger.info(f"Collected {len(uids)} UIDs total")
        return uids

    def fetch_metadata(
        self,
        uids: list[str],
    ) -> pd.DataFrame:
        """Fetch run-level metadata for a list of SRA UIDs.

        Uses raw XML parsing (ElementTree) since Biopython's Entrez.read()
        cannot handle SRA XML (no DTD).
        """
        all_rows: list[dict] = []

        for start in range(0, len(uids), BATCH_SIZE):
            batch = uids[start : start + BATCH_SIZE]
            handle = Entrez.efetch(
                db=SRA_DB,
                id=",".join(batch),
                rettype="full",
                retmode="xml",
            )
            raw_xml = handle.read()
            handle.close()

            try:
                root = ET.parse(BytesIO(raw_xml)).getroot()
            except ET.ParseError as e:
                logger.error(f"XML parse error at batch {start}: {e}")
                time.sleep(self.delay)
                continue

            for pkg in root.findall("EXPERIMENT_PACKAGE"):
                rows = self._parse_experiment_package(pkg)
                all_rows.extend(rows)

            logger.info(
                f"  parsed metadata {start}–{start + len(batch)} "
                f"({len(all_rows)} runs so far)"
            )
            time.sleep(self.delay)

        df = pd.DataFrame(all_rows)
        logger.info(f"Metadata collected for {len(df)} runs")
        return df

    def search_and_collect(
        self,
        query: str = DEFAULT_SRA_QUERY,
        max_results: int = 0,
        output: Optional[Path] = None,
    ) -> pd.DataFrame:
        """Convenience: search + fetch metadata, optionally save to TSV."""
        uids = self.search(query=query, max_results=max_results)
        df = self.fetch_metadata(uids)
        if output:
            output = Path(output)
            output.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output, sep="\t", index=False)
            logger.info(f"Metadata saved to {output}")
        return df

    # ------------------------------------------------------------------
    # Internal parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_experiment_package(pkg: ET.Element) -> list[dict]:
        """Extract fields from one EXPERIMENT_PACKAGE XML element.

        Returns a list of dicts (one per run in the package).
        """
        rows = []
        try:
            exp = pkg.find("EXPERIMENT")
            if exp is None:
                return rows

            # Experiment-level
            srx = _attr(exp, "accession")
            title = _text(exp, "TITLE")

            # Platform / instrument
            platform = ""
            instrument = ""
            plat_elem = exp.find("PLATFORM")
            if plat_elem is not None and len(plat_elem) > 0:
                plat_child = plat_elem[0]
                platform = plat_child.tag
                instrument = _text(plat_child, "INSTRUMENT_MODEL")

            # Library info
            lib = exp.find("DESIGN/LIBRARY_DESCRIPTOR")
            strategy = ""
            source = ""
            layout = ""
            if lib is not None:
                strategy = _text(lib, "LIBRARY_STRATEGY")
                source = _text(lib, "LIBRARY_SOURCE")
                layout_elem = lib.find("LIBRARY_LAYOUT")
                if layout_elem is not None and len(layout_elem) > 0:
                    layout = layout_elem[0].tag

            # Sample-level
            sample = pkg.find("SAMPLE")
            srs = _attr(sample, "accession") if sample is not None else ""
            organism = ""
            strain = ""
            if sample is not None:
                organism = _text(sample, "SAMPLE_NAME/SCIENTIFIC_NAME")
                for sa in sample.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
                    tag = _text(sa, "TAG").lower()
                    if tag == "strain":
                        strain = _text(sa, "VALUE")

            # Study-level
            study = pkg.find("STUDY")
            srp = _attr(study, "accession") if study is not None else ""
            bioproject = ""
            if study is not None:
                for ext_id in study.findall("IDENTIFIERS/EXTERNAL_ID"):
                    if ext_id.get("namespace") == "BioProject":
                        bioproject = ext_id.text or ""

            # Run-level
            for run in pkg.findall("RUN_SET/RUN"):
                srr = _attr(run, "accession")
                total_bases = _attr(run, "total_bases")
                total_spots = _attr(run, "total_spots")

                rows.append({
                    "srr_accession": srr,
                    "srx_accession": srx,
                    "srs_accession": srs,
                    "srp_accession": srp,
                    "bioproject": bioproject,
                    "organism": organism,
                    "strain": strain,
                    "platform": platform,
                    "instrument": instrument,
                    "total_bases": total_bases,
                    "total_spots": total_spots,
                    "library_strategy": strategy,
                    "library_source": source,
                    "library_layout": layout,
                    "title": title,
                })

        except Exception as e:
            logger.warning(f"Failed to parse experiment package: {e}")

        return rows
