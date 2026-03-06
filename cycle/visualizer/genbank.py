"""
IS Element GenBank Export Module (Cycle Pipeline)

Generates annotated GenBank (.gbk) files for IS elements with flanking regions,
ORFs, noncoding regions, and guide alignment hit positions. Output files are
loadable in SnapGene, Benchling, and other sequence viewers.

Adapted from RNA_guide_editor_finder/modules/is_element_genbank.py for the
Cycle JSON format (*_is_records_guide.json).
"""

import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor
from glob import glob
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .visualizer import _deduplicate_alignments

logger = logging.getLogger(__name__)


class ISElementGenBank:
    """Generate annotated GenBank files for IS elements."""

    def build_record(
        self,
        element: Dict,
        alignments: List[Dict],
        circle_info: Optional[Dict] = None,
        partial_circles: Optional[List[Dict]] = None,
    ) -> SeqRecord:
        """Build a BioPython SeqRecord for one IS element with annotated features.

        Coordinate system matches the visualizer: upstream flanking first,
        then IS body, then downstream flanking.

        Args:
            element: dict from *_is_records_guide.json.
            alignments: list of guide_hits for this element.
            circle_info: optional dict with circle_evidence fields.
            partial_circles: optional list of partial circle calls.

        Returns:
            Bio.SeqRecord.SeqRecord with features.
        """
        is_id = element.get("is_id", "unknown")
        is_data = element.get("is_element", {})
        is_seq = is_data.get("sequence", "")
        is_length = is_data.get("length", len(is_seq))

        up_seq = element.get("flanking_upstream", {}).get("sequence", "")
        down_seq = element.get("flanking_downstream", {}).get("sequence", "")
        up_len = len(up_seq)
        down_len = len(down_seq)

        full_seq = up_seq + is_seq + down_seq

        orf_ann = element.get("orf_annotation", {})
        orfs = orf_ann.get("orfs", [])
        nc_regions = orf_ann.get("noncoding_regions", [])

        features = []

        # 1. Flanking regions — misc_feature
        if up_len > 0:
            features.append(SeqFeature(
                FeatureLocation(0, up_len, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": ["upstream_flanking"],
                    "note": [f"{up_len}bp upstream flanking region"],
                },
            ))
        if down_len > 0:
            features.append(SeqFeature(
                FeatureLocation(up_len + is_length, up_len + is_length + down_len, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": ["downstream_flanking"],
                    "note": [f"{down_len}bp downstream flanking region"],
                },
            ))

        # 2. IS element boundary — misc_feature
        is_note = f"{is_length}bp IS element"
        sample = element.get("sample_id", "")
        if sample:
            is_note += f" sample={sample}"
        features.append(SeqFeature(
            FeatureLocation(up_len, up_len + is_length, strand=0),
            type="misc_feature",
            qualifiers={
                "label": [is_id],
                "note": [is_note],
            },
        ))

        # 2b. Circular IS annotation — misc_feature (if circle_info provided)
        if circle_info:
            th = circle_info.get("n_tail_head_reads", 0)
            gh = circle_info.get("n_genome_head_reads", 0)
            tg = circle_info.get("n_tail_genome_reads", 0)
            circle_note = (
                f"n_tail_head_reads={th}; "
                f"n_genome_head_reads={gh}; "
                f"n_tail_genome_reads={tg}"
            )
            features.append(SeqFeature(
                FeatureLocation(up_len, up_len + is_length, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": ["circular_IS"],
                    "note": [circle_note],
                },
            ))

        # 3. ORFs — CDS (with domain annotations if present)
        for i, orf in enumerate(orfs):
            orf_start = up_len + orf["start"] - 1  # 1-based to 0-based + offset
            orf_end = up_len + orf["end"]
            strand = +1 if orf.get("strand", "+") == "+" else -1
            length_nt = orf.get("length_nt", orf_end - orf_start)
            length_aa = length_nt // 3
            protein_seq = orf.get("protein_sequence", "")
            domains = orf.get("domains", [])
            orf_id = f"orf_{i+1}"

            qualifiers = {
                "label": [orf_id],
                "note": [f"{length_aa}aa"],
            }
            if protein_seq:
                qualifiers["translation"] = [protein_seq]
            if domains:
                qualifiers["note"] = [f"{length_aa}aa; domains: {','.join(domains)}"]
                qualifiers["product"] = [",".join(domains)]

            features.append(SeqFeature(
                FeatureLocation(orf_start, orf_end, strand=strand),
                type="CDS",
                qualifiers=qualifiers,
            ))

        # 4. Noncoding regions — misc_feature
        for nc in nc_regions:
            nc_start = up_len + nc["start"] - 1  # 1-based to 0-based + offset
            nc_end = up_len + nc["end"]
            nc_type = nc.get("type", "noncoding")
            features.append(SeqFeature(
                FeatureLocation(nc_start, nc_end, strand=0),
                type="misc_feature",
                qualifiers={
                    "label": [nc_type],
                    "note": [f"{nc_end - nc_start}bp noncoding region"],
                },
            ))

        # 5. Guide hits — misc_binding for both noncoding and flanking sides
        deduped = _deduplicate_alignments(alignments)
        for aln in deduped:
            nc_start_1based = aln.get("noncoding_start", 0)
            pos_in_nc = aln.get("pos_in_noncoding", 0)
            pos_in_flank = aln.get("pos_in_flanking", 0)
            aln_len = aln.get("length", 0)

            flanking = aln.get("flanking_source", "?")
            orientation = aln.get("orientation", "forward")
            mismatches = aln.get("mismatches", 0)

            # Build note
            note_parts = [f"orientation={orientation}", f"flanking_source={flanking}"]
            if mismatches:
                note_parts.append(f"mismatches={mismatches}")
            note = "; ".join(note_parts)

            label = f"aln_{flanking}_{aln_len}bp"

            # Noncoding-side hit
            nc_hit_start = up_len + (nc_start_1based - 1) + pos_in_nc
            nc_hit_end = nc_hit_start + aln_len
            features.append(SeqFeature(
                FeatureLocation(nc_hit_start, nc_hit_end, strand=0),
                type="misc_binding",
                qualifiers={
                    "label": [label + "_nc"],
                    "note": [note],
                },
            ))

            # Flanking-side hit
            if flanking == "upstream" and up_len > 0:
                flank_hit_start = pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            elif flanking == "downstream" and down_len > 0:
                flank_hit_start = up_len + is_length + pos_in_flank
                flank_hit_end = flank_hit_start + aln_len
            else:
                continue

            features.append(SeqFeature(
                FeatureLocation(flank_hit_start, flank_hit_end, strand=0),
                type="misc_binding",
                qualifiers={
                    "label": [label + "_flank"],
                    "note": [note],
                },
            ))

        # 6. Partial circle regions
        if partial_circles:
            for pc in partial_circles:
                pc_start = up_len + pc["circle_start"]
                pc_end = up_len + pc["circle_end"]
                n_reads = pc.get("n_supporting_reads", 0)
                frac = pc.get("circle_fraction", 0)
                pc_size = pc.get("circle_size", pc_end - pc_start)
                features.append(SeqFeature(
                    FeatureLocation(pc_start, pc_end, strand=0),
                    type="misc_feature",
                    qualifiers={
                        "label": [f"partial_circle_{pc_size}bp"],
                        "note": [
                            f"partial circle: {pc_size}bp ({frac:.0%} of IS), "
                            f"{n_reads} supporting reads"
                        ],
                    },
                ))

        # Build record
        description = f"{is_id} {is_length}bp IS element"

        record = SeqRecord(
            Seq(full_seq),
            id=is_id,
            name=is_id[:16],
            description=description,
            features=features,
            annotations={
                "molecule_type": "DNA",
                "topology": "linear",
            },
        )

        return record

    def save_genbank(
        self,
        element: Dict,
        alignments: List[Dict],
        output_path: str,
        circle_info: Optional[Dict] = None,
        partial_circles: Optional[List[Dict]] = None,
    ):
        """Generate and save a GenBank file for one IS element.

        Args:
            element: dict from *_is_records_guide.json.
            alignments: guide_hits for this element.
            output_path: path to save the .gbk file.
            circle_info: optional dict with circle_evidence fields.
            partial_circles: optional list of partial circle calls.
        """
        record = self.build_record(element, alignments, circle_info, partial_circles=partial_circles)
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        SeqIO.write(record, output_path, "genbank")

    def export_sample(
        self,
        json_path: str,
        output_dir: str,
    ):
        """Process one *_is_records_guide.json file, writing a GBK per record.

        Args:
            json_path: path to the guide JSON file.
            output_dir: directory to write GBK files into (a sample subdirectory is created).
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
            output_path = os.path.join(sample_dir, f"{is_id}.gbk")
            try:
                self.save_genbank(rec, alignments, output_path, circle_info=circle_info)
            except Exception:
                logger.exception("Failed to export GenBank for %s/%s", sample_id, is_id)

        logger.info("Exported %d GenBank records for %s", len(records), sample_id)

    def export_batch(
        self,
        formatter_dir: str,
        output_dir: str,
        parallel: int = 1,
    ):
        """Process all *_is_records_guide.json files under formatter_dir.

        Args:
            formatter_dir: directory containing sample subdirectories with guide JSONs.
            output_dir: directory to write GBK files.
            parallel: number of parallel workers.
        """
        json_files = sorted(glob(os.path.join(formatter_dir, "*", "*_is_records_guide.json")))
        if not json_files:
            logger.warning("No *_is_records_guide.json files found under %s", formatter_dir)
            return

        logger.info("Found %d guide JSON files to export", len(json_files))

        if parallel <= 1:
            for jf in json_files:
                self.export_sample(jf, output_dir)
        else:
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = [
                    pool.submit(_export_sample_worker, jf, output_dir)
                    for jf in json_files
                ]
                for fut in futures:
                    fut.result()


def _export_sample_worker(json_path, output_dir):
    """Top-level function for ProcessPoolExecutor (must be picklable)."""
    ISElementGenBank().export_sample(json_path, output_dir)
