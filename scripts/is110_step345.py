import csv
import json
import glob
import logging
import os
import sys
from collections import defaultdict

csv.field_size_limit(10_000_000)
sys.path.insert(0, '/home/kuangh/tools/Cycle')
from cycle.is110_filter import IS110Filter

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = "/groups/rubin/projects/kuang/out/IS_cycle"
OUT = f"{BASE}/is110_all_006_026"

# Reload all elements from Sniffles tables
logger.info("Loading Sniffles tables...")
all_elements = {}
for b in range(6, 27):
    snif_dir = f"{BASE}/batch_{b:03d}/sniffles_output"
    if not os.path.isdir(snif_dir):
        continue
    for sample_dir in sorted(os.listdir(snif_dir)):
        table_file = os.path.join(snif_dir, sample_dir, f"{sample_dir}.table.txt")
        if not os.path.exists(table_file):
            continue
        with open(table_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                uuid = row["UUID"]
                consensus = row.get("Consensus", "")
                if not consensus or len(consensus) < 500:
                    continue
                all_elements[uuid] = {
                    "sample_id": sample_dir,
                    "batch": f"batch_{b:03d}",
                    "chrom": row["Chrom"],
                    "start": int(row["Start"]),
                    "end": int(row["End"]),
                    "consensus": consensus,
                    "length": len(consensus),
                }
logger.info("Total IS elements: %d", len(all_elements))

# Parse Prodigal proteins
logger.info("Parsing Prodigal output...")
uuid_orfs = defaultdict(list)
n_proteins = 0
with open(f"{OUT}/all_proteins.faa") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header and seq_lines:
                seq = "".join(seq_lines).rstrip("*")
                uuid_orfs[header["uuid"]].append({
                    "prot_id": header["prot_id"],
                    "start": header["start"],
                    "end": header["end"],
                    "strand": header["strand"],
                    "sequence": seq,
                })
                n_proteins += 1
            parts = line[1:].strip().split(" # ")
            if len(parts) < 4:
                header = None
                seq_lines = []
                continue
            prot_id = parts[0]
            uuid = "_".join(prot_id.rsplit("_", 1)[:-1])
            strand = "+" if parts[3] == "1" else "-"
            header = {
                "prot_id": prot_id, "uuid": uuid,
                "start": int(parts[1]), "end": int(parts[2]), "strand": strand,
            }
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header and seq_lines:
        seq = "".join(seq_lines).rstrip("*")
        uuid_orfs[header["uuid"]].append({
            "prot_id": header["prot_id"], "start": header["start"],
            "end": header["end"], "strand": header["strand"], "sequence": seq,
        })
        n_proteins += 1
logger.info("Prodigal: %d proteins from %d elements", n_proteins, len(uuid_orfs))

# Write HMM-ready FASTA
hmm_fasta = f"{OUT}/all_proteins_hmm.faa"
with open(hmm_fasta, "w") as f:
    for uuid, orfs in uuid_orfs.items():
        for orf in orfs:
            hdr = f"{uuid}__{orf['start']}_{orf['end']}_{orf['strand']}"
            f.write(f">{hdr}\n{orf['sequence']}\n")
logger.info("Wrote HMM FASTA: %d proteins", n_proteins)

# HMM search
filt = IS110Filter()
is110_ids, trans_domains = filt.run_hmmsearch(hmm_fasta, OUT, cpus=4)
logger.info("IS110 elements: %d", len(is110_ids))

# Load circle evidence
logger.info("Loading circle evidence...")
circle_evidence = {}
for b in range(6, 27):
    circle_dirs = glob.glob(f"{BASE}/batch_{b:03d}/chunks/chunk_*/circle_output")
    for cdir in circle_dirs:
        if not os.path.isdir(cdir):
            continue
        for sample_dir in os.listdir(cdir):
            summary = os.path.join(cdir, sample_dir, f"{sample_dir}_circle_summary.tsv")
            if not os.path.exists(summary):
                continue
            with open(summary) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    uuid = row.get("is_uuid", "")
                    if uuid:
                        circle_evidence[uuid] = {
                            "n_tail_head_reads": int(row.get("n_tail_head_reads", 0)),
                            "n_genome_head_reads": int(row.get("n_genome_head_reads", 0)),
                            "n_tail_genome_reads": int(row.get("n_tail_genome_reads", 0)),
                            "n_total_mapped": int(row.get("n_total_mapped", 0)),
                        }
logger.info("Circle evidence for %d elements", len(circle_evidence))

pc_by_uuid = {}
for b in range(6, 27):
    pc_dir = f"{BASE}/batch_{b:03d}/partial_circle_output"
    for pcf in glob.glob(f"{pc_dir}/**/*_partial_circle_summary.json", recursive=True):
        with open(pcf) as f:
            data = json.load(f)
        if isinstance(data, list):
            for call in data:
                is_id = call.get("is_id", "")
                if is_id:
                    pc_by_uuid.setdefault(is_id, []).append(call)
logger.info("Partial circle data for %d elements", len(pc_by_uuid))

# Build records and split
with_circle = []
without_circle = []

for uuid in is110_ids:
    if uuid not in all_elements:
        continue
    elem = all_elements[uuid]
    ce = circle_evidence.get(uuid, {})
    n_th = ce.get("n_tail_head_reads", 0)
    has_partial = uuid in pc_by_uuid

    orfs = uuid_orfs.get(uuid, [])
    domains_info = trans_domains.get(uuid, {})
    for orf in orfs:
        orf_domains = []
        for domain_name in ("DEDD", "Tnp20"):
            for _, s, e, _ in domains_info.get(domain_name, []):
                if s == orf["start"] and e == orf["end"]:
                    orf_domains.append(domain_name)
        if orf_domains:
            orf["domains"] = sorted(set(orf_domains))

    rec = {
        "is_id": uuid,
        "sample_id": elem["sample_id"],
        "is_element": {
            "sequence": elem["consensus"],
            "length": elem["length"],
            "chrom": elem["chrom"],
            "start": elem["start"],
            "end": elem["end"],
        },
        "orf_annotation": {
            "num_orfs": len(orfs),
            "orfs": [{
                "start": o["start"], "end": o["end"], "strand": o["strand"],
                "length_nt": o["end"] - o["start"] + 1,
                "protein_sequence": o["sequence"],
                **({"domains": o["domains"]} if "domains" in o else {}),
            } for o in orfs],
        },
        "circle_evidence": ce if ce else {
            "n_tail_head_reads": 0, "n_genome_head_reads": 0,
            "n_tail_genome_reads": 0, "n_total_mapped": 0,
        },
        "_batch": elem["batch"],
    }

    if n_th > 0 or has_partial:
        rec["_circle_type"] = "full" if n_th > 0 else "partial"
        with_circle.append(rec)
    else:
        rec["_circle_type"] = "none"
        without_circle.append(rec)

logger.info("With circle: %d, Without: %d, Total: %d",
            len(with_circle), len(without_circle), len(with_circle) + len(without_circle))

filt.export_results(with_circle, without_circle, OUT)

# Visualize
from cycle.visualizer.visualizer import ISElementVisualizer
from cycle.visualizer.genbank import ISElementGenBank

vis = ISElementVisualizer()
gb = ISElementGenBank()

for label, records in [("with_circle_evidence", with_circle),
                       ("without_circle_evidence", without_circle)]:
    sub_dir = os.path.join(OUT, label)
    n_done = n_err = 0
    for rec in records:
        tag = f"{rec['sample_id']}__{rec['is_id'][:8]}"
        png_path = os.path.join(sub_dir, f"{tag}.png")
        gbk_path = os.path.join(sub_dir, f"{tag}.gbk")
        alignments = rec.get("guide_hits", [])
        circle_info = rec.get("circle_evidence")
        partial_circles = pc_by_uuid.get(rec["is_id"])
        try:
            vis.save_element_png(rec, alignments, png_path, dpi=150,
                                 circle_info=circle_info, partial_circles=partial_circles)
            gb.save_genbank(rec, alignments, gbk_path, circle_info=circle_info)
            n_done += 1
        except Exception as e:
            logger.error("Failed %s: %s", tag, e)
            n_err += 1
        if n_done % 100 == 0 and n_done > 0:
            logger.info("  %s: %d/%d done", label, n_done, len(records))
    logger.info("%s: %d PNG+GBK, %d errors", label, n_done, n_err)

logger.info("Done.")
