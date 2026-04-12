"""Plot IS110 split insertion: [IS110] vs [IS110 part1][unknown][IS110 part2]

Uses pygenomeviz ribbon style. Auto-detects whether to flip the genome track
so ribbons don't cross.
"""
import csv
import os
import sys

sys.path.insert(0, '/home/kuangh/tools/Cycle')

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)


def load_domains(hmm_dir, is110_ids):
    """Load DEDD/Tnp20 domain annotations from HMM tblout files.

    Returns dict: is110_id -> {(start, end): [domain_names]}
    """
    from collections import defaultdict
    from cycle.is110_filter import IS110Filter

    filt = IS110Filter()
    # coord_domains[uuid][(start,end)] = [domain_names]
    coord_domains = defaultdict(lambda: defaultdict(list))

    for name in ("DEDD", "Tnp20"):
        tbl = os.path.join(hmm_dir, f"{name}_hits.tbl")
        if not os.path.exists(tbl):
            continue
        hits = filt.parse_tblout(tbl)
        for prot_id, (trans_id, start, end, strand) in hits.items():
            if trans_id in is110_ids:
                coord_domains[trans_id][(start, end)].append(name)

    return dict(coord_domains)


def load_splits(tsv_path, max_insertion=10000):
    """Load split hits with insertion < max_insertion."""
    splits = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if int(row["insertion_length"]) <= max_insertion:
                splits.append(row)
    return splits


def get_orfs_batch(protein_fasta, is110_ids):
    """Get ORF coordinates from Prodigal protein FASTA for multiple IS110s."""
    from collections import defaultdict
    uuid_orfs = defaultdict(list)
    need = set(is110_ids)
    with open(protein_fasta) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split(" # ")
                if len(parts) < 4:
                    continue
                prot_id = parts[0]
                uuid = "_".join(prot_id.rsplit("_", 1)[:-1])
                if uuid in need:
                    uuid_orfs[uuid].append({
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "strand": 1 if parts[3] == "1" else -1,
                    })
    return uuid_orfs


def should_flip(split_info):
    """Detect if genome track should be flipped so ribbons don't cross.

    If part1 on query is AFTER part2 on query, but part1 is BEFORE part2 on ref,
    the ribbons will cross. Flipping the genome track fixes this.
    """
    qs1 = int(split_info["query_start_1"])
    qs2 = int(split_info["query_start_2"])
    rs1 = int(split_info["ref_start_1"])
    rs2 = int(split_info["ref_start_2"])

    # Check if order is reversed between query and reference
    query_order = qs1 < qs2  # True if part1 comes first on query
    ref_order = rs1 < rs2    # True if part1 comes first on ref

    return query_order != ref_order


def plot_split(split_info, orfs, output_prefix, orf_domains=None):
    """Create comparison: IS110 circular form vs genome split form.

    Args:
        orf_domains: dict mapping (start, end) -> [domain_names] for this IS110.
    """
    from pygenomeviz import GenomeViz

    if orf_domains is None:
        orf_domains = {}

    qlen = int(split_info["query_length"])
    qs1, qe1 = int(split_info["query_start_1"]), int(split_info["query_end_1"])
    qs2, qe2 = int(split_info["query_start_2"]), int(split_info["query_end_2"])
    rs1, re1 = int(split_info["ref_start_1"]), int(split_info["ref_end_1"])
    rs2, re2 = int(split_info["ref_start_2"]), int(split_info["ref_end_2"])
    ins_len = int(split_info["insertion_length"])
    pid1 = float(split_info["pident_1"])
    pid2 = float(split_info["pident_2"])
    contig = split_info["contig"]

    flip = should_flip(split_info)

    # Target region
    pad = 200
    tstart = rs1 - pad
    tend = re2 + pad
    tlen = tend - tstart

    # ── Build plot ──
    gv = GenomeViz(
        fig_width=16,
        fig_track_height=0.8,
        track_align_type="center",
        feature_track_ratio=0.5,
        link_track_ratio=1.2,
    )

    # Track 1: IS110 element (circular form)
    is110_name = f"IS110 ({qlen} bp) — circular form"
    track1 = gv.add_feature_track(is110_name, qlen)

    track1.add_feature(0, qlen, strand=1, plotstyle="bigarrow",
                       fc="#E8F5E9", ec="black", lw=1.0,
                       label=f"IS110 ({qlen}bp)", text_kws=dict(size=8))

    # Add ORFs with domain annotation
    DOMAIN_COLORS = {
        "DEDD+Tnp20": "#EE6677",  # red — both domains
        "DEDD": "#CC79A7",         # pink — DEDD only
        "Tnp20": "#E69F00",        # orange — Tnp20 only
        "other": "#4477AA",        # blue — no IS110 domains
    }

    for orf in orfs:
        domains = orf_domains.get((orf["start"], orf["end"]), [])
        if "DEDD" in domains and "Tnp20" in domains:
            color = DOMAIN_COLORS["DEDD+Tnp20"]
            label = "DEDD+Tnp20"
        elif "DEDD" in domains:
            color = DOMAIN_COLORS["DEDD"]
            label = "DEDD"
        elif "Tnp20" in domains:
            color = DOMAIN_COLORS["Tnp20"]
            label = "Tnp20"
        else:
            color = DOMAIN_COLORS["other"]
            label = ""

        track1.add_feature(orf["start"], orf["end"], strand=orf["strand"],
                           plotstyle="arrow", fc=color, ec="black", lw=0.5,
                           label=label, text_kws=dict(size=7, color="black"))

    # Mark parts on query
    track1.add_feature(qs1, qe1, strand=1, plotstyle="box",
                       fc="#BBDEFB", ec="#1565C0", lw=1.0, alpha=0.3,
                       label=f"part1 ({qs1}-{qe1})", text_kws=dict(size=7))
    track1.add_feature(qs2, qe2, strand=1, plotstyle="box",
                       fc="#C8E6C9", ec="#2E7D32", lw=1.0, alpha=0.3,
                       label=f"part2 ({qs2}-{qe2})", text_kws=dict(size=7))

    # Track 2: Genome region
    if flip:
        genome_name = f"{contig} ({tend:,}→{tstart:,}) [rev comp]"
    else:
        genome_name = f"{contig} ({tstart:,}→{tend:,})"
    track2 = gv.add_feature_track(genome_name, tlen)

    def ref_to_track(pos):
        """Convert reference position to track coordinate, handling flip."""
        t = pos - tstart
        if flip:
            return tlen - t
        return t

    # Part 1 on genome
    t_rs1 = ref_to_track(rs1)
    t_re1 = ref_to_track(re1)
    p1_left, p1_right = min(t_rs1, t_re1), max(t_rs1, t_re1)

    # Insertion
    t_ins_start = ref_to_track(re1)
    t_ins_end = ref_to_track(rs2)
    ins_left, ins_right = min(t_ins_start, t_ins_end), max(t_ins_start, t_ins_end)

    # Part 2 on genome
    t_rs2 = ref_to_track(rs2)
    t_re2 = ref_to_track(re2)
    p2_left, p2_right = min(t_rs2, t_re2), max(t_rs2, t_re2)

    # Determine display order after potential flip
    if p1_left < p2_left:
        # Part 1 is left, part 2 is right
        track2.add_feature(p1_left, p1_right, strand=1, plotstyle="bigarrow",
                           fc="#BBDEFB", ec="#1565C0", lw=1.0,
                           label=f"IS110 part1 ({re1-rs1}bp)", text_kws=dict(size=7))
        track2.add_feature(ins_left, ins_right, strand=1, plotstyle="box",
                           fc="#FFCDD2", ec="#C62828", lw=1.5,
                           label=f"unknown ({ins_len}bp)", text_kws=dict(size=8, color="red"))
        track2.add_feature(p2_left, p2_right, strand=1, plotstyle="bigarrow",
                           fc="#C8E6C9", ec="#2E7D32", lw=1.0,
                           label=f"IS110 part2 ({re2-rs2}bp)", text_kws=dict(size=7))
    else:
        # Part 2 is left, part 1 is right (after flip)
        track2.add_feature(p2_left, p2_right, strand=1, plotstyle="bigarrow",
                           fc="#C8E6C9", ec="#2E7D32", lw=1.0,
                           label=f"IS110 part2 ({re2-rs2}bp)", text_kws=dict(size=7))
        track2.add_feature(ins_left, ins_right, strand=1, plotstyle="box",
                           fc="#FFCDD2", ec="#C62828", lw=1.5,
                           label=f"unknown ({ins_len}bp)", text_kws=dict(size=8, color="red"))
        track2.add_feature(p1_left, p1_right, strand=1, plotstyle="bigarrow",
                           fc="#BBDEFB", ec="#1565C0", lw=1.0,
                           label=f"IS110 part1 ({re1-rs1}bp)", text_kws=dict(size=7))

    # ── Links ──
    link1_target = (p1_left, p1_right) if not flip else (min(t_rs1, t_re1), max(t_rs1, t_re1))
    link2_target = (p2_left, p2_right) if not flip else (min(t_rs2, t_re2), max(t_rs2, t_re2))

    gv.add_link(
        (is110_name, qs1, qe1),
        (genome_name, p1_left, p1_right),
        color="blue", v=pid1, vmin=80, vmax=100, curve=True, alpha=0.5,
    )
    gv.add_link(
        (is110_name, qs2, qe2),
        (genome_name, p2_left, p2_right),
        color="blue", v=pid2, vmin=80, vmax=100, curve=True, alpha=0.5,
    )

    gv.set_colorbar(
        colors=["blue"], vmin=80, vmax=100,
        bar_label="Identity (%)",
        bar_height=0.15, bar_width=0.01,
        bar_labelsize=10, tick_labelsize=8,
    )

    gv.savefig(f"{output_prefix}.png", dpi=200)


def main():
    BASE = "/groups/rubin/projects/kuang/out/IS_cycle"
    SPLIT_TSV = f"{BASE}/is110_split_gtdb/is110_split_hits.tsv"
    PROTEIN_FASTA = f"{BASE}/is110_all_006_026/all_proteins.faa"
    OUT_DIR = f"{BASE}/is110_split_gtdb/examples"
    os.makedirs(OUT_DIR, exist_ok=True)

    # Load splits with insertion < 10000bp
    splits = load_splits(SPLIT_TSV, max_insertion=10000)
    logger.info("Loaded %d splits with insertion < 10kb", len(splits))

    # Limit to 1000
    splits = splits[:1000]
    logger.info("Generating PNGs for %d splits", len(splits))

    # Batch load all ORFs and domains
    is110_ids = list(set(s["is110_id"] for s in splits))
    logger.info("Loading ORFs for %d unique IS110s...", len(is110_ids))
    uuid_orfs = get_orfs_batch(PROTEIN_FASTA, is110_ids)
    logger.info("Loaded ORFs for %d IS110s", len(uuid_orfs))

    HMM_DIR = f"{BASE}/is110_all_006_026"
    logger.info("Loading domain annotations...")
    all_domains = load_domains(HMM_DIR, set(is110_ids))
    logger.info("Loaded domains for %d IS110s", len(all_domains))

    n_done = 0
    n_err = 0
    for split_info in splits:
        is110_id = split_info["is110_id"]
        contig = split_info["contig"].replace("|", "_").replace("/", "_")
        tag = f"{is110_id[:8]}_{contig[:20]}_{split_info['insertion_length']}bp"
        output_prefix = os.path.join(OUT_DIR, tag)

        orfs = uuid_orfs.get(is110_id, [])
        orf_domains = all_domains.get(is110_id, {})
        try:
            plot_split(split_info, orfs, output_prefix, orf_domains=orf_domains)
            n_done += 1
        except Exception as e:
            logger.error("Failed %s: %s", tag, e)
            n_err += 1

        if n_done % 50 == 0 and n_done > 0:
            logger.info("  %d/%d done", n_done, len(splits))

    logger.info("Done: %d PNGs, %d errors", n_done, n_err)


if __name__ == "__main__":
    main()
