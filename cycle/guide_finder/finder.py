"""Guide finder — search for short alignments between flanking and noncoding regions."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

from .alignment import ShortAlignmentFinder

logger = logging.getLogger(__name__)


class GuideFinder:
    """Find potential RNA-guide alignments between flanking and noncoding regions.

    For each IS element with ORF annotations, aligns flanking sequences
    against internal noncoding regions.  Alignments suggest the element
    encodes a guide that recognises its target site.

    Uses tiered mismatch rule:
      - Exact (0 mm): min_length (default 9bp)
      - With mismatches: min_length_for_mismatch (default 12bp)

    Parameters
    ----------
    min_length : int
        Minimum alignment length for exact matches (default 9).
    max_mismatches : int
        Maximum mismatches allowed for longer hits (default 1).
    min_length_for_mismatch : int
        Minimum alignment length when mismatches > 0 (default 12).
    check_revcomp : bool
        Also search reverse-complement orientations (default True).
    """

    def __init__(
        self,
        min_length: int = 9,
        max_mismatches: int = 1,
        min_length_for_mismatch: int = 12,
        check_revcomp: bool = True,
    ) -> None:
        self.aligner = ShortAlignmentFinder(
            min_length=min_length,
            max_mismatches=max_mismatches,
            min_length_for_mismatch=min_length_for_mismatch,
            check_forward=True,
            check_revcomp=check_revcomp,
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def find_guides_sample(self, json_path: str | Path) -> Optional[Path]:
        """Run guide search on a single annotated JSON file.

        Parameters
        ----------
        json_path : path to ``*_is_records_annotated.json``.

        Returns
        -------
        Path to ``*_is_records_guide.json``, or None on failure.
        """
        json_path = Path(json_path)
        if not json_path.exists():
            logger.warning(f"File not found: {json_path}")
            return None

        with open(json_path) as fh:
            records = json.load(fh)

        if not records:
            logger.warning(f"No records in {json_path.name}")
            return None

        total_hits = 0
        for rec in records:
            hits = self._search_record(rec)
            rec["guide_hits"] = hits
            rec["guide_summary"] = self._summarise(hits)
            total_hits += len(hits)

        out_path = json_path.with_name(
            json_path.name.replace(
                "_is_records_annotated.json", "_is_records_guide.json"
            )
        )
        with open(out_path, "w") as fh:
            json.dump(records, fh, indent=2)

        sample_id = json_path.stem.replace("_is_records_annotated", "")
        logger.info(
            f"{sample_id}: {total_hits} guide hits across {len(records)} records -> {out_path.name}"
        )
        return out_path

    def find_guides_batch(
        self, formatter_dir: str | Path, parallel: int = 1,
        sample_ids: set[str] | None = None,
    ) -> dict[str, Optional[Path]]:
        """Run guide search on all annotated JSONs under *formatter_dir*.

        Glob pattern: ``*/*_is_records_annotated.json`` (one level deep).
        If *sample_ids* is given, only process those sample directories.
        """
        formatter_dir = Path(formatter_dir)
        json_files = sorted(
            formatter_dir.glob("*/*_is_records_annotated.json")
        )

        if sample_ids is not None:
            json_files = [jf for jf in json_files if jf.parent.name in sample_ids]

        if not json_files:
            logger.warning(
                f"No *_is_records_annotated.json files in {formatter_dir}"
            )
            return {}

        logger.info(f"Found {len(json_files)} annotated JSONs in {formatter_dir}")

        results: dict[str, Optional[Path]] = {}

        if parallel <= 1:
            for jf in json_files:
                try:
                    results[jf.name] = self.find_guides_sample(jf)
                except Exception as exc:
                    logger.error(f"Failed for {jf.name}: {exc}")
                    results[jf.name] = None
        else:
            from concurrent.futures import ProcessPoolExecutor, as_completed

            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for jf in json_files:
                    fut = pool.submit(
                        _guide_worker,
                        jf,
                        self.aligner.min_length,
                        self.aligner.max_mismatches,
                        self.aligner.min_length_for_mismatch,
                        self.aligner.check_revcomp,
                    )
                    futures[fut] = jf.name

                for fut in as_completed(futures):
                    name = futures[fut]
                    try:
                        results[name] = fut.result()
                    except Exception as exc:
                        logger.error(f"Worker failed for {name}: {exc}")
                        results[name] = None

        succeeded = sum(1 for v in results.values() if v is not None)
        logger.info(f"Batch done: {succeeded}/{len(results)} samples processed")
        return results

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _search_record(self, rec: dict) -> list[dict]:
        """Search one IS record for guide alignments."""
        orf_ann = rec.get("orf_annotation")
        if not orf_ann:
            return []

        noncoding = orf_ann.get("noncoding_regions")
        if not noncoding:
            return []

        is_seq = (rec.get("is_element") or {}).get("sequence", "")
        if not is_seq:
            return []

        flanks = []
        up = (rec.get("flanking_upstream") or {}).get("sequence", "")
        if up:
            flanks.append(("upstream", up))
        down = (rec.get("flanking_downstream") or {}).get("sequence", "")
        if down:
            flanks.append(("downstream", down))

        if not flanks:
            return []

        hits: list[dict] = []

        for nc_idx, nc in enumerate(noncoding):
            # Noncoding coordinates are 1-based from Prodigal
            nc_start = nc["start"]
            nc_end = nc["end"]
            nc_seq = is_seq[nc_start - 1 : nc_end]

            if not nc_seq:
                continue

            for flank_source, flank_seq in flanks:
                alignments = self.aligner.find_alignments_between(
                    flank_seq, nc_seq
                )
                for aln in alignments:
                    hits.append(
                        {
                            "flanking_source": flank_source,
                            "noncoding_region_index": nc_idx,
                            "noncoding_region_type": nc["type"],
                            "noncoding_start": nc_start,
                            "noncoding_end": nc_end,
                            "pos_in_flanking": aln["pos1"],
                            "pos_in_noncoding": aln["pos2"],
                            "length": aln["length"],
                            "orientation": aln["orientation"],
                            "mismatches": aln["mismatches"],
                            "mismatch_positions": aln["mismatch_positions"],
                            "seq_flanking": aln["seq1"],
                            "seq_noncoding": aln["seq2"],
                        }
                    )

        return hits

    @staticmethod
    def _summarise(hits: list[dict]) -> dict:
        if not hits:
            return {"n_hits": 0, "best_length": 0, "has_revcomp_hit": False}
        return {
            "n_hits": len(hits),
            "best_length": max(h["length"] for h in hits),
            "has_revcomp_hit": any(
                h["orientation"] == "reverse_complement" for h in hits
            ),
        }


# ------------------------------------------------------------------
# Worker for parallel execution
# ------------------------------------------------------------------


def _guide_worker(
    json_path: Path,
    min_length: int,
    max_mismatches: int,
    min_length_for_mismatch: int,
    check_revcomp: bool,
) -> Optional[Path]:
    """Standalone worker for ProcessPoolExecutor."""
    gf = GuideFinder(
        min_length=min_length,
        max_mismatches=max_mismatches,
        min_length_for_mismatch=min_length_for_mismatch,
        check_revcomp=check_revcomp,
    )
    return gf.find_guides_sample(json_path)
