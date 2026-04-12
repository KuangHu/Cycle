"""Collect all 2276 IS110 candidates, merge with circle evidence from formatter output."""
import csv
import json
import glob
import logging
import os
import sys
from collections import defaultdict

csv.field_size_limit(10_000_000)
sys.path.insert(0, '/home/kuangh/tools/Cycle')

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = "/groups/rubin/projects/kuang/out/IS_cycle"
HMM_DIR = f"{BASE}/is110_all_006_026"
OUT = f"{BASE}/is110_all_006_026/collected"
os.makedirs(OUT, exist_ok=True)

# Load IS110 IDs
with open(f"{HMM_DIR}/is110_ids.json") as f:
    is110_ids = set(json.load(f))
logger.info("IS110 IDs: %d", len(is110_ids))

# Step 1: Collect full records from formatter output (have flanking, guide, circle evidence)
logger.info("Step 1: Collecting from formatter guide JSONs...")
formatter_records = {}  # is_id -> record
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
                rec["_source"] = "formatter"
                formatter_records[is_id] = rec
logger.info("From formatter (with TH): %d", len(formatter_records))

# Step 2: Load Prodigal ORFs for remaining IS110s
logger.info("Step 2: Parsing Prodigal ORFs...")
uuid_orfs = defaultdict(list)
with open(f"{HMM_DIR}/all_proteins.faa") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header and seq_lines:
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
            strand = "+" if parts[3] == "1" else "-"
            header = {"uuid": uuid, "start": int(parts[1]), "end": int(parts[2]), "strand": strand}
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header and seq_lines:
        seq = "".join(seq_lines).rstrip("*")
        uuid_orfs[header["uuid"]].append({
            "start": header["start"], "end": header["end"],
            "strand": header["strand"], "protein_sequence": seq,
            "length_nt": header["end"] - header["start"] + 1,
        })

# Step 3: Load Sniffles tables for consensus sequences of non-formatter IS110s
logger.info("Step 3: Loading Sniffles consensus for remaining IS110s...")
missing_ids = is110_ids - set(formatter_records.keys())
logger.info("Need Sniffles data for %d IS110s", len(missing_ids))

sniffles_records = {}
for b in range(6, 29):
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
                if uuid not in missing_ids:
                    continue
                consensus = row.get("Consensus", "")
                orfs = uuid_orfs.get(uuid, [])
                sniffles_records[uuid] = {
                    "is_id": uuid,
                    "sample_id": sample_dir,
                    "is_element": {
                        "sequence": consensus,
                        "length": len(consensus),
                        "chrom": row["Chrom"],
                        "start": int(row["Start"]),
                        "end": int(row["End"]),
                    },
                    "orf_annotation": {
                        "num_orfs": len(orfs),
                        "orfs": orfs,
                    },
                    "circle_evidence": {
                        "n_tail_head_reads": 0,
                        "n_genome_head_reads": 0,
                        "n_tail_genome_reads": 0,
                        "n_total_mapped": 0,
                    },
                    "guide_hits": [],
                    "guide_summary": {"n_hits": 0, "best_length": 0, "has_revcomp_hit": False},
                    "_batch": f"batch_{b:03d}",
                    "_source": "sniffles",
                }
logger.info("From Sniffles (no TH): %d", len(sniffles_records))

# Step 4: Load circle evidence from circle_summary TSVs (some might have reads but no TH)
logger.info("Step 4: Loading circle evidence...")
for b in range(6, 29):
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
                    if uuid in sniffles_records:
                        sniffles_records[uuid]["circle_evidence"] = {
                            "n_tail_head_reads": int(row.get("n_tail_head_reads", 0)),
                            "n_genome_head_reads": int(row.get("n_genome_head_reads", 0)),
                            "n_tail_genome_reads": int(row.get("n_tail_genome_reads", 0)),
                            "n_total_mapped": int(row.get("n_total_mapped", 0)),
                        }

# Step 5: Annotate domains on all records
logger.info("Step 5: Annotating domains...")
from cycle.is110_filter import IS110Filter
filt = IS110Filter()
protein_hits = {}
for name in ("DEDD", "Tnp20"):
    tbl = f"{HMM_DIR}/{name}_hits.tbl"
    protein_hits[name] = filt.parse_tblout(tbl)
_, trans_domains = filt._identify_is110(protein_hits)

for rec in list(formatter_records.values()) + list(sniffles_records.values()):
    filt._annotate_orf_domains(rec, trans_domains)

# Step 6: Split and export
logger.info("Step 6: Splitting and exporting...")
with_circle = []
without_circle = []

all_records = list(formatter_records.values()) + list(sniffles_records.values())
for rec in all_records:
    ce = rec.get("circle_evidence", {})
    n_th = ce.get("n_tail_head_reads", 0)
    if n_th > 0:
        rec["_circle_type"] = "full"
        with_circle.append(rec)
    else:
        rec["_circle_type"] = "none"
        without_circle.append(rec)

logger.info("With circle evidence: %d", len(with_circle))
logger.info("Without circle evidence: %d", len(without_circle))
logger.info("Total: %d", len(all_records))

# Export
for label, records in [("with_circle_evidence", with_circle),
                       ("without_circle_evidence", without_circle)]:
    sub_dir = os.path.join(OUT, label)
    os.makedirs(sub_dir, exist_ok=True)

    with open(os.path.join(sub_dir, "is110_records.json"), "w") as f:
        json.dump(records, f, indent=2)

    fields = [
        "is_id", "sample_id", "batch", "source", "is_length", "n_orfs", "circle_type",
        "n_tail_head_reads", "n_genome_head_reads", "n_tail_genome_reads",
        "n_guide_hits", "best_guide_length",
    ]
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
                "source": rec.get("_source", ""),
                "is_length": rec.get("is_element", {}).get("length", 0),
                "n_orfs": (rec.get("orf_annotation") or {}).get("num_orfs", 0),
                "circle_type": rec.get("_circle_type", ""),
                "n_tail_head_reads": ce.get("n_tail_head_reads", 0),
                "n_genome_head_reads": ce.get("n_genome_head_reads", 0),
                "n_tail_genome_reads": ce.get("n_tail_genome_reads", 0),
                "n_guide_hits": gs.get("n_hits", 0),
                "best_guide_length": gs.get("best_length", 0),
            })
    logger.info("Wrote %s: %d records", label, len(records))

logger.info("Done. Output: %s", OUT)
