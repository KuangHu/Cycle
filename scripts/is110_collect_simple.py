"""Collect all IS110 candidates, split into with/without circle evidence.

- With circle: full records from formatter guide JSONs (612)
- Without circle: IS110 IDs not in formatter (need consensus from Sniffles for GBK)
  Only loads Sniffles data for the ~1600 missing IS110s, not all 774K.
"""
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
HMM_DIR = f"{BASE}/is110_all_006_026"
OUT = f"{BASE}/is110_all_006_026/collected"

# Load IS110 IDs
with open(f"{HMM_DIR}/is110_ids.json") as f:
    is110_ids = set(json.load(f))
logger.info("IS110 IDs: %d", len(is110_ids))

# Load domain annotations
filt = IS110Filter()
protein_hits = {}
for name in ("DEDD", "Tnp20"):
    protein_hits[name] = filt.parse_tblout(f"{HMM_DIR}/{name}_hits.tbl")
_, trans_domains = filt._identify_is110(protein_hits)

# Step 1: Collect from formatter guide JSONs (with circle evidence)
logger.info("Collecting from formatter guide JSONs...")
with_circle = []
with_circle_ids = set()
for b in range(6, 29):
    fmt_dir = f"{BASE}/batch_{b:03d}/is_formatter_output"
    if not os.path.isdir(fmt_dir):
        continue
    for jf in sorted(glob.glob(f"{fmt_dir}/*/*_is_records_guide.json")):
        with open(jf) as f:
            records = json.load(f)
        for rec in records:
            is_id = rec.get("is_id", "")
            if is_id in is110_ids:
                rec["_batch"] = f"batch_{b:03d}"
                rec["_circle_type"] = "full"
                filt._annotate_orf_domains(rec, trans_domains)
                with_circle.append(rec)
                with_circle_ids.add(is_id)
logger.info("With circle evidence: %d", len(with_circle))

# Step 2: For remaining IS110s, build index of which batch/sample they're in
missing_ids = is110_ids - with_circle_ids
logger.info("Without circle evidence: %d to find", len(missing_ids))

# Build sample lookup from Sniffles dirs — only read tables that might have our IDs
# First build batch/sample index from sniffles_output dir structure
without_circle = []
found = 0
for b in range(6, 29):
    snif_dir = f"{BASE}/batch_{b:03d}/sniffles_output"
    if not os.path.isdir(snif_dir):
        continue
    for sample_dir in os.listdir(snif_dir):
        table_file = os.path.join(snif_dir, sample_dir, f"{sample_dir}.table.txt")
        if not os.path.exists(table_file):
            continue
        with open(table_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                uuid = row["UUID"]
                if uuid not in missing_ids:
                    continue
                consensus = row.get("Consensus", "")
                # Get Prodigal ORFs from HMM fasta headers
                rec = {
                    "is_id": uuid,
                    "sample_id": sample_dir,
                    "is_element": {
                        "sequence": consensus,
                        "length": len(consensus),
                        "chrom": row["Chrom"],
                        "start": int(row["Start"]),
                        "end": int(row["End"]),
                    },
                    "circle_evidence": {
                        "n_tail_head_reads": 0,
                        "n_genome_head_reads": 0,
                        "n_tail_genome_reads": 0,
                        "n_total_mapped": 0,
                    },
                    "_batch": f"batch_{b:03d}",
                    "_circle_type": "none",
                }
                filt._annotate_orf_domains(rec, trans_domains)
                without_circle.append(rec)
                found += 1
        if found >= len(missing_ids):
            break
    if found >= len(missing_ids):
        break

logger.info("Found %d/%d without circle evidence", len(without_circle), len(missing_ids))

# Parse Prodigal ORFs for without_circle records
logger.info("Adding Prodigal ORFs to without-circle records...")
uuid_orfs = defaultdict(list)
need_orfs = {r["is_id"] for r in without_circle}
with open(f"{HMM_DIR}/all_proteins.faa") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header and header["uuid"] in need_orfs and seq_lines:
                seq = "".join(seq_lines).rstrip("*")
                uuid_orfs[header["uuid"]].append({
                    "start": header["start"], "end": header["end"],
                    "strand": header["strand"], "protein_sequence": seq,
                    "length_nt": header["end"] - header["start"] + 1,
                })
            parts = line[1:].strip().split(" # ")
            if len(parts) < 4:
                header = None
                seq_lines = []
                continue
            prot_id = parts[0]
            uuid = "_".join(prot_id.rsplit("_", 1)[:-1])
            if uuid not in need_orfs:
                header = None
                seq_lines = []
                continue
            strand = "+" if parts[3] == "1" else "-"
            header = {"uuid": uuid, "start": int(parts[1]), "end": int(parts[2]), "strand": strand}
            seq_lines = []
        else:
            if header:
                seq_lines.append(line.strip())
    if header and header["uuid"] in need_orfs and seq_lines:
        seq = "".join(seq_lines).rstrip("*")
        uuid_orfs[header["uuid"]].append({
            "start": header["start"], "end": header["end"],
            "strand": header["strand"], "protein_sequence": seq,
            "length_nt": header["end"] - header["start"] + 1,
        })

for rec in without_circle:
    orfs = uuid_orfs.get(rec["is_id"], [])
    rec["orf_annotation"] = {"num_orfs": len(orfs), "orfs": orfs}
    filt._annotate_orf_domains(rec, trans_domains)

# Export
logger.info("Exporting...")
from cycle.visualizer.genbank import ISElementGenBank
gb = ISElementGenBank()

for label, records in [("with_circle_evidence", with_circle),
                       ("without_circle_evidence", without_circle)]:
    sub_dir = os.path.join(OUT, label)
    os.makedirs(sub_dir, exist_ok=True)

    with open(os.path.join(sub_dir, "is110_records.json"), "w") as f:
        json.dump(records, f, indent=2)

    # Summary TSV
    fields = ["is_id", "sample_id", "batch", "is_length", "n_orfs", "circle_type",
              "n_tail_head_reads", "n_genome_head_reads", "n_tail_genome_reads",
              "n_guide_hits", "best_guide_length"]
    with open(os.path.join(sub_dir, "is110_summary.tsv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rec in records:
            ce = rec.get("circle_evidence", {})
            gs = rec.get("guide_summary", {})
            writer.writerow({
                "is_id": rec.get("is_id", ""),
                "sample_id": rec.get("sample_id", ""),
                "batch": rec.get("_batch", ""),
                "is_length": rec.get("is_element", {}).get("length", 0),
                "n_orfs": (rec.get("orf_annotation") or {}).get("num_orfs", 0),
                "circle_type": rec.get("_circle_type", ""),
                "n_tail_head_reads": ce.get("n_tail_head_reads", 0),
                "n_genome_head_reads": ce.get("n_genome_head_reads", 0),
                "n_tail_genome_reads": ce.get("n_tail_genome_reads", 0),
                "n_guide_hits": gs.get("n_hits", 0),
                "best_guide_length": gs.get("best_length", 0),
            })

    # GenBank files
    n_done = n_err = 0
    for rec in records:
        tag = f"{rec['sample_id']}__{rec['is_id'][:8]}"
        gbk_path = os.path.join(sub_dir, f"{tag}.gbk")
        try:
            gb.save_genbank(rec, rec.get("guide_hits", []), gbk_path,
                            circle_info=rec.get("circle_evidence"))
            n_done += 1
        except Exception as e:
            logger.error("Failed %s: %s", tag, e)
            n_err += 1
    logger.info("%s: %d GBK, %d errors", label, n_done, n_err)

logger.info("Done. Output: %s", OUT)
