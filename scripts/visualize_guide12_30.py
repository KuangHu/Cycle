#!/usr/bin/env python
"""
Generate PNG + GBK for all IS elements in the top 7 clusters
whose best guide alignment length is between 12 and 30 (inclusive).
"""

import json
import logging
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from cycle.visualizer import ISElementVisualizer, ISElementGenBank

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

FORMATTER_DIR = "/groups/rubin/projects/kuang/out/IS_cycle/batch_000/is_formatter_output"
CLUSTER_TSV = "/groups/rubin/projects/kuang/out/IS_cycle/system_clustering_batch000/system_clusters_summary.tsv"
OUTPUT_DIR = "/groups/rubin/projects/kuang/out/IS_cycle/visualizations_guide12_30"
TARGET_CLUSTERS = {0, 1, 3, 4, 5, 10, 13}
MIN_GL = 12
MAX_GL = 30


def load_target_elements():
    """Read cluster TSV and return set of (sample_id, is_id) for targets, plus cluster mapping."""
    targets = {}  # is_id -> cluster_id
    samples_needed = set()

    with open(CLUSTER_TSV) as f:
        header = f.readline().strip().split("\t")
        ci = {col: i for i, col in enumerate(header)}

        for line in f:
            fields = line.strip().split("\t")
            cluster_id = int(fields[ci["cluster_id"]])
            if cluster_id not in TARGET_CLUSTERS:
                continue

            best_gl = int(fields[ci["best_guide_length"]])
            if best_gl < MIN_GL or best_gl > MAX_GL:
                continue

            is_id = fields[ci["is_id"]]
            sample_id = fields[ci["sample_id"]]
            targets[is_id] = cluster_id
            samples_needed.add(sample_id)

    return targets, samples_needed


def main():
    logger.info("Loading cluster assignments...")
    targets, samples_needed = load_target_elements()
    logger.info("Found %d target elements across %d samples", len(targets), len(samples_needed))

    # Load guide JSONs for needed samples, index by is_id
    records_by_id = {}
    for sample_id in sorted(samples_needed):
        json_path = os.path.join(FORMATTER_DIR, sample_id, f"{sample_id}_is_records_guide.json")
        if not os.path.exists(json_path):
            logger.warning("Missing guide JSON for %s", sample_id)
            continue
        with open(json_path) as f:
            records = json.load(f)
        for rec in records:
            rid = rec.get("is_id")
            if rid in targets:
                records_by_id[rid] = rec

    logger.info("Loaded %d / %d target records from guide JSONs", len(records_by_id), len(targets))

    vis = ISElementVisualizer()
    gb = ISElementGenBank()

    done = 0
    failed = 0
    for is_id, rec in sorted(records_by_id.items(), key=lambda kv: targets[kv[0]]):
        cluster_id = targets[is_id]
        sample_id = rec.get("sample_id", "unknown")
        alignments = rec.get("guide_hits", [])
        circle_info = rec.get("circle_evidence")

        out_dir = os.path.join(OUTPUT_DIR, f"cluster_{cluster_id}", sample_id)
        os.makedirs(out_dir, exist_ok=True)

        png_path = os.path.join(out_dir, f"{is_id}.png")
        gbk_path = os.path.join(out_dir, f"{is_id}.gbk")

        try:
            vis.save_element_png(rec, alignments, png_path, circle_info=circle_info)
            gb.save_genbank(rec, alignments, gbk_path, circle_info=circle_info)
            done += 1
        except Exception:
            logger.exception("Failed: cluster=%d sample=%s is_id=%s", cluster_id, sample_id, is_id)
            failed += 1

        if done % 20 == 0:
            logger.info("Progress: %d done, %d failed", done, failed)

    logger.info("Complete: %d done, %d failed out of %d total", done, failed, done + failed)


if __name__ == "__main__":
    main()
