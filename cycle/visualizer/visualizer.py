"""
IS Element Visualizer Module (Cycle Pipeline)

Generates per-element PNG diagrams using dna_features_viewer showing IS element
structure: flanking regions, ORFs (directional arrows), noncoding regions, and
guide alignment hit positions.

Adapted from RNA_guide_editor_finder/modules/is_element_visualizer.py for the
Cycle JSON format (*_is_records_guide.json).
"""

import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

logger = logging.getLogger(__name__)

# Color palette for IS families
IS_FAMILY_COLORS = {
    "IS110": "#e41a1c",
    "IS3": "#377eb8",
    "IS5": "#4daf4a",
    "IS21": "#984ea3",
    "IS630": "#ff7f00",
    "IS1182": "#a65628",
    "IS481": "#f781bf",
    "IS66": "#999999",
    "IS256": "#66c2a5",
    "ISL3": "#fc8d62",
    "IS30": "#8da0cb",
    "IS4": "#e78ac3",
    "IS1": "#a6d854",
    "IS1595": "#ffd92f",
    "IS91": "#e5c494",
    "IS701": "#b3b3b3",
    "IS1380": "#1b9e77",
    "IS982": "#d95f02",
    "IS6": "#7570b3",
    "IS1634": "#e7298a",
    "ISNCY": "#66a61e",
    "unclassified": "#cccccc",
}
DEFAULT_ORF_COLOR = "#bbbbbb"
NC_COLOR = "#e0e0e0"
FLANKING_UP_COLOR = "#aed6f1"     # light blue for upstream flanking
FLANKING_DOWN_COLOR = "#f9e79f"   # light yellow for downstream flanking
ALIGNMENT_UPSTREAM_COLOR = "#d62728"    # red for upstream flanking hits
ALIGNMENT_DOWNSTREAM_COLOR = "#1f77b4"  # blue for downstream flanking hits


def _deduplicate_alignments(alignments: List[Dict]) -> List[Dict]:
    """Remove duplicate alignments based on position, length, source, and orientation."""
    seen = set()
    deduped = []
    for aln in alignments:
        key = (
            aln.get("flanking_source"),
            aln.get("noncoding_start"),
            aln.get("pos_in_noncoding"),
            aln.get("length"),
            aln.get("orientation"),
        )
        if key not in seen:
            seen.add(key)
            deduped.append(aln)
    return deduped


class ISElementVisualizer:
    """Generate diagrams of IS elements with flanking regions, ORFs, noncoding regions, and guide hits."""

    def visualize_element(
        self,
        element: Dict,
        alignments: List[Dict],
    ) -> plt.Figure:
        """Build a dna_features_viewer diagram for one IS element.

        The diagram spans from -upstream_length to is_length + downstream_length,
        showing flanking regions on either side of the IS element body.

        Args:
            element: dict from *_is_records_guide.json.
            alignments: list of guide_hits for this element.

        Returns:
            matplotlib Figure.
        """
        is_length = element.get("is_element", {}).get("length", 0)
        orf_ann = element.get("orf_annotation", {})
        orfs = orf_ann.get("orfs", [])
        nc_regions = orf_ann.get("noncoding_regions", [])

        # Flanking region lengths
        up_len = element.get("flanking_upstream", {}).get("length", 0)
        down_len = element.get("flanking_downstream", {}).get("length", 0)

        # Coordinate system: flanking_upstream occupies [-up_len, 0),
        # IS element occupies [0, is_length), downstream occupies [is_length, is_length+down_len)
        total_length = up_len + is_length + down_len
        offset = up_len  # shift all IS-internal coords by this amount

        features = []

        # 0. Flanking regions
        if up_len > 0:
            features.append(GraphicFeature(
                start=0,
                end=up_len,
                strand=0,
                color=FLANKING_UP_COLOR,
                label=f"upstream ({up_len}bp)",
                linewidth=0.5,
            ))
        if down_len > 0:
            features.append(GraphicFeature(
                start=offset + is_length,
                end=offset + is_length + down_len,
                strand=0,
                color=FLANKING_DOWN_COLOR,
                label=f"downstream ({down_len}bp)",
                linewidth=0.5,
            ))

        # 1. Noncoding regions — light gray rectangles
        for nc in nc_regions:
            nc_start = offset + nc["start"] - 1  # 1-based to 0-based + offset
            nc_end = offset + nc["end"]
            nc_type = nc.get("type", "")
            label = nc_type.replace("_", " ") if nc_type else "NC"
            nc_len = nc_end - nc_start
            features.append(GraphicFeature(
                start=nc_start,
                end=nc_end,
                strand=0,
                color=NC_COLOR,
                label=label if nc_len >= 30 else None,
                linewidth=0.5,
            ))

        # 2. ORFs — colored directional arrows
        for i, orf in enumerate(orfs):
            orf_start = offset + orf["start"] - 1  # 1-based to 0-based + offset
            orf_end = offset + orf["end"]
            strand = +1 if orf.get("strand", "+") == "+" else -1
            length_nt = orf.get("length_nt", orf_end - orf_start)
            length_aa = length_nt // 3

            label = f"orf_{i+1} ({length_aa}aa)"

            features.append(GraphicFeature(
                start=orf_start,
                end=orf_end,
                strand=strand,
                color=DEFAULT_ORF_COLOR,
                label=label,
                linewidth=1,
            ))

        # 3. Guide hits — deduplicated, colored markers on both flanking and noncoding sides
        deduped = _deduplicate_alignments(alignments)
        alignment_pairs = []  # collect (nc_midpoint, flank_midpoint, color) for connecting lines
        for aln in deduped:
            nc_start_1based = aln.get("noncoding_start", 0)
            pos_in_nc = aln.get("pos_in_noncoding", 0)
            pos_in_flank = aln.get("pos_in_flanking", 0)
            aln_len = aln.get("length", 0)

            flanking = aln.get("flanking_source", "?")
            color = ALIGNMENT_UPSTREAM_COLOR if flanking == "upstream" else ALIGNMENT_DOWNSTREAM_COLOR
            orientation = aln.get("orientation", "forward")
            arrow = "\u2191" if flanking == "upstream" else "\u2193"
            ori_tag = "" if orientation == "forward" else " rc"
            label = f"{arrow}{aln_len}bp{ori_tag}"

            # Noncoding-side position (0-based diagram coords)
            nc_hit_start = offset + (nc_start_1based - 1) + pos_in_nc
            nc_hit_end = nc_hit_start + aln_len

            features.append(GraphicFeature(
                start=nc_hit_start,
                end=nc_hit_end,
                strand=0,
                color=color,
                label=label,
                linewidth=1.5,
                linecolor=color,
            ))

            # Flanking-side position (0-based diagram coords)
            if flanking == "upstream" and up_len > 0:
                flank_hit_start = pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            elif flanking == "downstream" and down_len > 0:
                flank_hit_start = offset + is_length + pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            else:
                continue  # no flanking region to draw in

            features.append(GraphicFeature(
                start=flank_hit_start,
                end=flank_hit_end,
                strand=0,
                color=color,
                label=None,
                linewidth=1.5,
                linecolor=color,
            ))

            # Save midpoints for connecting lines
            nc_mid = (nc_hit_start + nc_hit_end) / 2
            flank_mid = (flank_hit_start + flank_hit_end) / 2
            alignment_pairs.append((nc_mid, flank_mid, color))

        record = GraphicRecord(sequence_length=total_length, features=features)
        # Cap figure width to avoid matplotlib pixel limit (65535 px max)
        fig_width = max(8, min(300, total_length / 150))
        ax, _ = record.plot(figure_width=fig_width)

        # Draw connecting lines between flanking-side and noncoding-side hits
        for nc_mid, flank_mid, color in alignment_pairs:
            ax.annotate(
                "", xy=(nc_mid, -0.4), xytext=(flank_mid, -0.4),
                arrowprops=dict(
                    arrowstyle="-",
                    color=color,
                    alpha=0.35,
                    linewidth=1.5,
                    connectionstyle="arc3,rad=0.3",
                ),
            )

        # Draw IS element boundary lines
        ax.axvline(x=offset, color="black", linewidth=1.5, linestyle="--", alpha=0.6)
        ax.axvline(x=offset + is_length, color="black", linewidth=1.5, linestyle="--", alpha=0.6)

        fig = ax.figure
        return fig

    def save_element_png(
        self,
        element: Dict,
        alignments: List[Dict],
        output_path: str,
        dpi: int = 150,
        figure_width: int = 12,
        circle_info: Optional[Dict] = None,
    ):
        """Generate and save a PNG diagram for one IS element.

        Args:
            element: dict from *_is_records_guide.json.
            alignments: guide_hits for this element.
            output_path: path to save the PNG.
            dpi: resolution.
            figure_width: figure width in inches.
            circle_info: optional dict with circle_evidence fields.
        """
        fig = self.visualize_element(element, alignments)

        # Build title
        is_id = element.get("is_id", "unknown")
        is_length = element.get("is_element", {}).get("length", 0)

        deduped_count = len(_deduplicate_alignments(alignments))
        title = f"{is_id}  ({is_length} bp)  \u2014  {deduped_count} guide hit(s)"

        if circle_info:
            th = circle_info.get("n_tail_head_reads", 0)
            gh = circle_info.get("n_genome_head_reads", 0)
            tg = circle_info.get("n_tail_genome_reads", 0)
            title += f"  |  circle: {th} TH, {gh} GH, {tg} TG reads"

        fig.axes[0].set_title(title, fontsize=10, pad=10)
        fig.set_size_inches(figure_width, fig.get_size_inches()[1])

        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    def visualize_sample(
        self,
        json_path: str,
        output_dir: str,
        dpi: int = 150,
        figure_width: int = 12,
    ):
        """Process one *_is_records_guide.json file, writing a PNG per record.

        Args:
            json_path: path to the guide JSON file.
            output_dir: directory to write PNGs into (a sample subdirectory is created).
            dpi: resolution.
            figure_width: figure width in inches.
        """
        with open(json_path) as f:
            records = json.load(f)

        if not records:
            logger.info("No records in %s, skipping", json_path)
            return

        sample_id = records[0].get("sample_id", "unknown")
        sample_dir = os.path.join(output_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)

        for rec in records:
            is_id = rec.get("is_id", "unknown")
            alignments = rec.get("guide_hits", [])
            circle_info = rec.get("circle_evidence")
            output_path = os.path.join(sample_dir, f"{is_id}.png")
            try:
                self.save_element_png(
                    rec, alignments, output_path,
                    dpi=dpi, figure_width=figure_width,
                    circle_info=circle_info,
                )
            except Exception:
                logger.exception("Failed to visualize %s/%s", sample_id, is_id)

        logger.info("Visualized %d records for %s", len(records), sample_id)

    def visualize_batch(
        self,
        formatter_dir: str,
        output_dir: str,
        parallel: int = 1,
        dpi: int = 150,
        figure_width: int = 12,
    ):
        """Process all *_is_records_guide.json files under formatter_dir.

        Args:
            formatter_dir: directory containing sample subdirectories with guide JSONs.
            output_dir: directory to write PNGs.
            parallel: number of parallel workers.
            dpi: resolution.
            figure_width: figure width in inches.
        """
        json_files = sorted(glob(os.path.join(formatter_dir, "*", "*_is_records_guide.json")))
        if not json_files:
            logger.warning("No *_is_records_guide.json files found under %s", formatter_dir)
            return

        logger.info("Found %d guide JSON files to visualize", len(json_files))

        if parallel <= 1:
            for jf in json_files:
                self.visualize_sample(jf, output_dir, dpi=dpi, figure_width=figure_width)
        else:
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = [
                    pool.submit(
                        _visualize_sample_worker, jf, output_dir, dpi, figure_width,
                    )
                    for jf in json_files
                ]
                for fut in futures:
                    fut.result()  # propagate exceptions


def _visualize_sample_worker(json_path, output_dir, dpi, figure_width):
    """Top-level function for ProcessPoolExecutor (must be picklable)."""
    ISElementVisualizer().visualize_sample(json_path, output_dir, dpi=dpi, figure_width=figure_width)
