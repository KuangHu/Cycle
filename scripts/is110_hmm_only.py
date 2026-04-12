"""Parse Prodigal output, write HMM-ready FASTA, run hmmsearch, report IS110 IDs."""
import json
import logging
import os
import sys
from collections import defaultdict

sys.path.insert(0, '/home/kuangh/tools/Cycle')
from cycle.is110_filter import IS110Filter

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

OUT = "/groups/rubin/projects/kuang/out/IS_cycle/is110_all_006_026"

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
                    "start": header["start"], "end": header["end"],
                    "strand": header["strand"], "sequence": seq,
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
            header = {"uuid": uuid, "start": int(parts[1]), "end": int(parts[2]), "strand": strand}
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header and seq_lines:
        seq = "".join(seq_lines).rstrip("*")
        uuid_orfs[header["uuid"]].append({
            "start": header["start"], "end": header["end"],
            "strand": header["strand"], "sequence": seq,
        })
        n_proteins += 1

logger.info("Parsed %d proteins from %d elements", n_proteins, len(uuid_orfs))

# Write HMM-ready FASTA
hmm_fasta = f"{OUT}/all_proteins_hmm.faa"
with open(hmm_fasta, "w") as f:
    for uuid, orfs in uuid_orfs.items():
        for orf in orfs:
            f.write(f">{uuid}__{orf['start']}_{orf['end']}_{orf['strand']}\n{orf['sequence']}\n")
logger.info("Wrote HMM FASTA: %s", hmm_fasta)

# Run hmmsearch
filt = IS110Filter()
is110_ids, trans_domains = filt.run_hmmsearch(hmm_fasta, OUT, cpus=44)

# Save IS110 IDs
with open(f"{OUT}/is110_ids.json", "w") as f:
    json.dump(list(is110_ids), f)
logger.info("Saved %d IS110 IDs to is110_ids.json", len(is110_ids))
