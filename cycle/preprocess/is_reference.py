"""Build a tldr-formatted IS element reference FASTA from ISfinder sequences."""

import logging
import re
import urllib.request
from pathlib import Path
from typing import Optional

from .config import DEFAULT_IS_REFERENCE_DIR, ISFINDER_FASTA_URL

logger = logging.getLogger(__name__)


class ISReferenceBuilder:
    """Download ISfinder sequences and reformat for tldr.

    ISfinder headers look like ``>ISname_ISgroup_ISfamily``.
    tldr expects ``>ISfamily:ISname``.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_IS_REFERENCE_DIR,
        url: str = ISFINDER_FASTA_URL,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.url = url

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(
        self,
        families: Optional[list[str]] = None,
        output_name: str = "is_reference.fa",
    ) -> Path:
        """Download, reformat, and (optionally) filter IS sequences.

        Args:
            families: If given, only keep sequences from these IS families
                      (e.g. ["IS256", "IS6", "IS3"]).  Case-insensitive.
            output_name: Output filename.

        Returns:
            Path to the output FASTA.
        """
        raw_path = self.output_dir / "ISfinder_raw.fna"
        out_path = self.output_dir / output_name

        # Download raw FASTA (skip if cached)
        self._download(raw_path)

        # Parse + reformat + optional filter
        records = self._parse_and_reformat(raw_path)

        if families:
            families_lower = {f.lower() for f in families}
            before = len(records)
            records = [
                r for r in records
                if r["family"].lower() in families_lower
            ]
            logger.info(
                f"Filtered to {len(records)}/{before} sequences "
                f"(families: {families})"
            )

        # Write output
        with open(out_path, "w") as fh:
            for rec in records:
                fh.write(f">{rec['family']}:{rec['name']}\n")
                fh.write(rec["seq"] + "\n")

        logger.info(
            f"IS reference written: {out_path} "
            f"({len(records)} sequences)"
        )
        return out_path

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _download(self, dest: Path) -> None:
        """Download raw ISfinder FASTA if not already cached."""
        if dest.exists():
            logger.info(f"Using cached ISfinder FASTA: {dest}")
            return

        logger.info(f"Downloading ISfinder sequences from {self.url}")
        urllib.request.urlretrieve(self.url, dest)
        logger.info(f"  → {dest}")

    @staticmethod
    def _parse_and_reformat(fasta: Path) -> list[dict]:
        """Parse ISfinder FASTA and extract name/family from headers.

        Expected header format: ``>ISname_ISgroup_ISfamily``
        (underscore-separated, at least 3 fields).

        Returns list of {name, group, family, seq}.
        """
        records: list[dict] = []
        current_header = ""
        seq_lines: list[str] = []

        def _flush():
            if not current_header:
                return
            info = _parse_header(current_header)
            if info:
                info["seq"] = "".join(seq_lines)
                records.append(info)

        with open(fasta) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    _flush()
                    current_header = line[1:].split()[0]  # first word only
                    seq_lines = []
                else:
                    seq_lines.append(line)
            _flush()

        logger.info(f"Parsed {len(records)} IS sequences")
        return records


def _parse_header(header: str) -> Optional[dict]:
    """Parse ``ISname_ISgroup_ISfamily`` into components.

    Handles edge cases:
      - Headers with fewer than 3 underscore-separated fields
      - Names containing underscores (take last two fields as group/family)
    """
    parts = header.split("_")
    if len(parts) >= 3:
        family = parts[-1]
        group = parts[-2]
        name = "_".join(parts[:-2])
        return {"name": name, "group": group, "family": family}

    # Fallback: try regex for IS-prefixed names
    m = re.match(r"^(IS\S+?)_(IS\S+?)_(IS\S+)$", header)
    if m:
        return {"name": m.group(1), "group": m.group(2), "family": m.group(3)}

    logger.warning(f"Cannot parse IS header: {header}")
    return None


# Make module-level logger available to the function
logger = logging.getLogger(__name__)
