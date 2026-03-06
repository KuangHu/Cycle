#!/usr/bin/env python3
"""Visualize nesting event pairs — PNG + GBK for host and core IS elements.

For each nesting event, outputs diagrams of both the host (longer) and core
(shorter) IS element, with:
  - Guide hits (flanking/noncoding alignments) as usual
  - Aligned blocks highlighted (green) on both host and core
  - Insertion regions highlighted (coral red) on the host

Example
-------
    python scripts/visualize_nesting_events.py \
        --events nesting_output/nesting_events.json \
        --records is110_circular_records.json \
        --output-dir nesting_output/visualizations
"""

import argparse
import json
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from cycle.visualizer.visualizer import ISElementVisualizer, _deduplicate_alignments
from cycle.visualizer.genbank import ISElementGenBank

logger = logging.getLogger(__name__)

# Colors
INSERTION_COLOR = "#ff6b6b"          # coral red for insertion regions
ALIGNED_COLOR = "#66bb6a"            # green for aligned blocks


# ---------------------------------------------------------------------------
# Extended-seq → diagram coordinate conversion
# ---------------------------------------------------------------------------

def _ext_to_diagram(ext_pos, flanking_pad, up_len):
    """Convert extended-sequence position to diagram coordinate.

    Extended seq = up_pad(flanking_pad) + IS + dn_pad(flanking_pad).
    Diagram       = upstream(up_len)   + IS + downstream(dn_len).
    Offset: diagram = ext_pos - flanking_pad + up_len
    """
    return ext_pos - flanking_pad + up_len


# ---------------------------------------------------------------------------
# PNG annotation helpers
# ---------------------------------------------------------------------------

def _annotate_blocks_png(fig, element, blocks, flanking_pad, coord_key):
    """Highlight aligned blocks on a PNG figure.

    Args:
        coord_key: 'host' or 'core' — which block coordinates to use.
    """
    ax = fig.axes[0]
    up_len = element.get("flanking_upstream", {}).get("length", 0)

    for i, blk in enumerate(blocks):
        start = _ext_to_diagram(blk[f"{coord_key}_start"], flanking_pad, up_len)
        end = _ext_to_diagram(blk[f"{coord_key}_end"], flanking_pad, up_len)

        ax.axvspan(start, end, alpha=0.15, color=ALIGNED_COLOR, zorder=0)
        ax.axvline(start, color=ALIGNED_COLOR, linewidth=1, linestyle=":",
                   alpha=0.6)
        ax.axvline(end, color=ALIGNED_COLOR, linewidth=1, linestyle=":",
                   alpha=0.6)

        mid = (start + end) / 2
        bp = end - start
        ax.text(mid, ax.get_ylim()[0] * 0.6, f"blk{i+1} {bp}bp",
                ha="center", va="top", fontsize=7, color="#2e7d32",
                fontstyle="italic")


def _annotate_insertions_png(fig, element, insertions, flanking_pad):
    """Highlight insertion regions on the host PNG figure."""
    ax = fig.axes[0]
    up_len = element.get("flanking_upstream", {}).get("length", 0)

    for ins in insertions:
        start = _ext_to_diagram(ins["host_start"], flanking_pad, up_len)
        end = _ext_to_diagram(ins["host_end"], flanking_pad, up_len)
        ins_size = ins["insertion_size"]

        ax.axvspan(start, end, alpha=0.25, color=INSERTION_COLOR, zorder=0)
        ax.axvline(start, color=INSERTION_COLOR, linewidth=1.5,
                   linestyle=":", alpha=0.8)
        ax.axvline(end, color=INSERTION_COLOR, linewidth=1.5,
                   linestyle=":", alpha=0.8)

        mid = (start + end) / 2
        ax.text(mid, ax.get_ylim()[1] * 0.85, f"ins {ins_size}bp",
                ha="center", va="bottom", fontsize=8, color=INSERTION_COLOR,
                fontweight="bold")


# ---------------------------------------------------------------------------
# GBK annotation helpers
# ---------------------------------------------------------------------------

def _annotate_blocks_gbk(record, element, blocks, flanking_pad, coord_key):
    """Add aligned block features to a BioPython SeqRecord."""
    up_len = len((element.get("flanking_upstream") or {}).get("sequence", ""))

    for i, blk in enumerate(blocks):
        start = _ext_to_diagram(blk[f"{coord_key}_start"], flanking_pad, up_len)
        end = _ext_to_diagram(blk[f"{coord_key}_end"], flanking_pad, up_len)
        start = max(0, start)
        end = min(len(record.seq), end)
        bp = end - start
        ident = blk.get("identity", 0)

        record.features.append(SeqFeature(
            FeatureLocation(start, end, strand=0),
            type="misc_feature",
            qualifiers={
                "label": [f"aligned_block_{i+1}"],
                "note": [f"{bp}bp aligned block; identity={ident:.1%}"],
            },
        ))


def _annotate_insertions_gbk(record, element, insertions, flanking_pad):
    """Add insertion region features to a BioPython SeqRecord."""
    up_len = len((element.get("flanking_upstream") or {}).get("sequence", ""))

    for i, ins in enumerate(insertions):
        start = _ext_to_diagram(ins["host_start"], flanking_pad, up_len)
        end = _ext_to_diagram(ins["host_end"], flanking_pad, up_len)
        start = max(0, start)
        end = min(len(record.seq), end)

        record.features.append(SeqFeature(
            FeatureLocation(start, end, strand=0),
            type="misc_feature",
            qualifiers={
                "label": [f"insertion_{i+1}"],
                "note": [f"{ins['insertion_size']}bp nested insertion"],
            },
        ))


# ---------------------------------------------------------------------------
# Per-event visualization
# ---------------------------------------------------------------------------

def visualize_event(event, records_by_id, vis, gb, output_dir, flanking_pad,
                    dpi=150, figure_width=14):
    """Generate PNG + GBK for one nesting event (host and core)."""
    host_id = event["host_is_id"]
    core_id = event["core_is_id"]
    host_rec = records_by_id.get(host_id)
    core_rec = records_by_id.get(core_id)

    if not host_rec or not core_rec:
        logger.warning("Record not found for %s or %s, skipping", host_id, core_id)
        return

    event_label = f"{host_id[:8]}_vs_{core_id[:8]}"
    event_dir = os.path.join(output_dir, event_label)
    os.makedirs(event_dir, exist_ok=True)

    insertions = event.get("insertions", [])
    blocks = event.get("aligned_blocks", [])

    # --- Host element ---
    host_alns = host_rec.get("guide_hits", [])
    host_circle = host_rec.get("circle_evidence")

    # PNG
    fig = vis.visualize_element(host_rec, host_alns)
    _annotate_blocks_png(fig, host_rec, blocks, flanking_pad, "host")
    _annotate_insertions_png(fig, host_rec, insertions, flanking_pad)

    deduped = len(_deduplicate_alignments(host_alns))
    title = (
        f"HOST: {host_id[:16]}...  ({event['host_is_length']}bp)  "
        f"-- {deduped} guide hit(s)  |  "
        f"{event['n_insertions']} insertion(s), {event['total_insertion_bp']}bp total"
    )
    if host_circle:
        th = host_circle.get("n_tail_head_reads", 0)
        title += f"  |  {th} TH reads"

    fig.axes[0].set_title(title, fontsize=9, pad=10)
    fig.set_size_inches(figure_width, fig.get_size_inches()[1])
    fig.savefig(os.path.join(event_dir, f"host_{host_id}.png"),
                dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    # GBK
    gbk_rec = gb.build_record(host_rec, host_alns, circle_info=host_circle)
    _annotate_blocks_gbk(gbk_rec, host_rec, blocks, flanking_pad, "host")
    _annotate_insertions_gbk(gbk_rec, host_rec, insertions, flanking_pad)
    SeqIO.write(gbk_rec, os.path.join(event_dir, f"host_{host_id}.gbk"), "genbank")

    # --- Core element ---
    core_alns = core_rec.get("guide_hits", [])
    core_circle = core_rec.get("circle_evidence")

    # PNG
    fig = vis.visualize_element(core_rec, core_alns)
    _annotate_blocks_png(fig, core_rec, blocks, flanking_pad, "core")

    deduped = len(_deduplicate_alignments(core_alns))
    title = (
        f"CORE: {core_id[:16]}...  ({event['core_is_length']}bp)  "
        f"-- {deduped} guide hit(s)"
    )
    if core_circle:
        th = core_circle.get("n_tail_head_reads", 0)
        title += f"  |  {th} TH reads"

    fig.axes[0].set_title(title, fontsize=9, pad=10)
    fig.set_size_inches(figure_width, fig.get_size_inches()[1])
    fig.savefig(os.path.join(event_dir, f"core_{core_id}.png"),
                dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    # GBK
    gbk_rec = gb.build_record(core_rec, core_alns, circle_info=core_circle)
    _annotate_blocks_gbk(gbk_rec, core_rec, blocks, flanking_pad, "core")
    SeqIO.write(gbk_rec, os.path.join(event_dir, f"core_{core_id}.gbk"), "genbank")

    logger.info("  %s: host=%dbp core=%dbp ins=%s blocks=%d",
                event_label,
                event["host_is_length"], event["core_is_length"],
                "+".join(str(i["insertion_size"]) for i in insertions),
                len(blocks))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--events", required=True,
        help="Path to nesting_events.json",
    )
    parser.add_argument(
        "--records", required=True,
        help="Path to is110 records JSON (is110_circular_records.json or similar)",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory for visualization output",
    )
    parser.add_argument(
        "--flanking-pad", type=int, default=80,
        help="Flanking pad used in nesting detection (default: 80)",
    )
    parser.add_argument(
        "--dpi", type=int, default=150,
        help="PNG resolution (default: 150)",
    )
    parser.add_argument(
        "--figure-width", type=int, default=14,
        help="Figure width in inches (default: 14)",
    )
    parser.add_argument(
        "--max-events", type=int, default=None,
        help="Maximum number of events to visualize (default: all)",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    events = json.load(open(args.events))
    records = json.load(open(args.records))
    records_by_id = {r["is_id"]: r for r in records}

    if args.max_events:
        events = events[:args.max_events]

    logger.info("Loaded %d events, %d records", len(events), len(records))

    vis = ISElementVisualizer()
    gb = ISElementGenBank()

    os.makedirs(args.output_dir, exist_ok=True)

    for i, event in enumerate(events):
        logger.info("Event %d/%d", i + 1, len(events))
        try:
            visualize_event(
                event, records_by_id, vis, gb, args.output_dir,
                flanking_pad=args.flanking_pad,
                dpi=args.dpi, figure_width=args.figure_width,
            )
        except Exception:
            logger.exception("Failed to visualize event %d", i + 1)

    logger.info("Done — %d events visualized -> %s", len(events), args.output_dir)


if __name__ == "__main__":
    main()
