# Part 1: Batch Pipeline

The batch pipeline runs on every batch (batch_000 through batch_094) via SLURM.
It takes raw SRA accessions and produces fully annotated IS element records with
circle evidence, ORF annotations, guide RNA hits, system clustering, novelty
scores, and partial circle calls.

## Overview

```
Step 1: Download + Resolve + Index     ─┐
Step 2: Align                           │  submit_all.sh
Step 3: Sniffles + Circle + Format      │  (one SLURM job per batch per step)
Step 4: ORF annotation + Guide finder   │
Step 5: Clustering + Novelty + Partial  ─┘
```

Each step can be chained with `--dep` so the next step waits for the previous.

## Quick Start

```bash
# Process batches 021-025, fully chained:
for s in 1 2 3 4 5; do
  dep=""; [ $s -gt 1 ] && dep="--dep $((s-1))"
  bash submit_all.sh --step $s $dep 21 25
done

# Or submit one step at a time:
bash submit_all.sh --step 1 21 25
# Wait, check logs, then:
bash submit_all.sh --step 2 --dep 1 21 25
# ...
```

## Steps in Detail

### Step 1: Download + Resolve + Index

**Scripts:** `run_pipeline.py --steps download resolve index`
**Resources:** 48 CPUs, 192G, 2 days
**Input:** `data/batches/batch_NNN.tsv` (SRA accessions + organism names)
**Output:**
- `batch_NNN/sra_downloads/` — FASTQ files (one per accession)
- `batch_NNN/reference_genomes/` — Reference genome FASTAs + minimap2 `.mmi` indices

**What it does:**
1. **download** — Fetches FASTQs via kingfisher (SRR/ERR/DRR prefixes)
2. **resolve** — Downloads one reference genome per organism via NCBI datasets CLI
3. **index** — Builds minimap2 `.mmi` index per reference genome

### Step 2: Align

**Scripts:** `run_pipeline.py --steps align`
**Resources:** 48 CPUs, 192G, 2 days
**Input:** FASTQs + reference indices from step 1
**Output:**
- `batch_NNN/alignments/{SRR}_{organism}.sorted.bam` — Sorted BAMs with `.bai` indices
- `batch_NNN/alignments/alignment_status.tsv` — Status log

**What it does:**
- Runs minimap2 (map-ont preset) + samtools sort per sample
- 44 threads for minimap2, 4G sort memory per thread

### Step 3: Sniffles + Circle + Format

**Scripts:** `prepare_batch_metadata.py` + `run_pipeline.py --steps sniffles circle format`
**Resources:** 48 CPUs, 192G, 2 days (6 parallel samples × 8 threads)
**Input:** Sorted BAMs from step 2
**Output:**
- `batch_NNN/sniffles_output/` — Sniffles2 SV calls (`.vcf`)
- `batch_NNN/circle_output/` — Concatemer bait junction detection results
- `batch_NNN/is_formatter_output/{SRR}/` — Per-sample IS element records:
  - `{SRR}_is_records.json` — Raw IS records with sequences + flanking
  - `{SRR}_is_records.tsv` — Flat summary

**What it does:**
1. **sniffles** — Calls structural variants (insertions) via Sniffles2
2. **circle** — Detects tail-head junction reads (full circle evidence) via concatemer bait
3. **format** — Extracts IS element sequences + flanking regions via local assembly

### Step 4: ORF Annotation + Guide Finder

**Scripts:** `run_orf_annotator.py` + `run_guide_finder.py`
**Resources:** 48 CPUs, 192G, 6 hours
**Input:** `batch_NNN/is_formatter_output/` from step 3
**Output:** Enriched JSON files in the same directory:
- `{SRR}_is_records_annotated.json` — Records with ORF annotations (Prodigal)
- `{SRR}_is_records_guide.json` — Records with guide RNA alignment hits

**What it does:**
1. **orf_annotator** — Runs Prodigal on each IS element, adds ORF coordinates + protein sequences + noncoding regions to each record
2. **guide_finder** — Aligns noncoding regions against flanking regions to find guide RNA candidates, adds guide_hits to each record

**Note:** The guide JSON (`_is_records_guide.json`) is the final per-sample record
and is the input for all downstream analysis (Part 2).

### Step 5: Clustering + Novelty + Partial Circle

**Scripts:** `run_system_clustering.py` + `run_novelty_annotator.py` + `run_partial_circle.py`
**Resources:** 48 CPUs, 192G, 12 hours
**Input:** Guide JSONs from step 4 + FASTQs from step 1
**Output:**
- `system_clustering_batch_NNN/` — System-level IS clustering (MMseqs2 + Louvain)
  - `system_clusters.json` — Cluster assignments
  - `all_proteins.faa` — All predicted proteins (used by IS110 filter)
- `novelty_batch_NNN/` — ISfinder novelty annotation
  - `cluster_novelty_summary.tsv` — Per-cluster novelty scores
- `batch_NNN/partial_circle_output/{SRR}/` — Per-sample partial circle calls
  - `{SRR}_partial_circle_summary.json` — Clustered back-jump breakpoints
  - `{SRR}_partial_circle_reads.tsv` — Raw split-read detections

**What it does:**
1. **clustering** — Clusters IS elements by shared protein families (MMseqs2 protein clustering → Louvain community detection → flanking/guide variant grouping)
2. **novelty** — BLASTs all IS consensus sequences against ISfinder database, assigns composite novelty score per cluster
3. **partial_circle** — Maps reads to single-copy IS references, detects split-read back-jumps indicating sub-element circularization

## Data Flow

```
batch_NNN.tsv
  │
  ├─ Step 1 ─→ sra_downloads/          (FASTQs)
  │             reference_genomes/       (refs + .mmi)
  │
  ├─ Step 2 ─→ alignments/             (sorted BAMs)
  │
  ├─ Step 3 ─→ sniffles_output/        (VCFs)
  │             circle_output/           (junction reads)
  │             is_formatter_output/     (_is_records.json per sample)
  │
  ├─ Step 4 ─→ is_formatter_output/     (_is_records_annotated.json)
  │             is_formatter_output/     (_is_records_guide.json)  ← final records
  │
  └─ Step 5 ─→ system_clustering_batch_NNN/  (clusters + all_proteins.faa)
               novelty_batch_NNN/             (novelty scores)
               partial_circle_output/         (back-jump calls per sample)
```

## Monitoring

```bash
# Check all running jobs
squeue -u $(whoami)

# Check specific batch
squeue -u $(whoami) -n batch_021_s3

# Live logs
tail -f ~/logs/batch_021_step3_*.err

# Check completion
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed
```

## Troubleshooting

- **Step 1 slow**: kingfisher retries on transient SRA failures. Check download_status.tsv.
- **Step 2 hangs**: Some large genomes take hours per sample. Check alignment_status.tsv.
- **Step 3 Sniffles error**: Ensure `--allow-overwrite` is set (already in pipeline). Check VCF has variants before skipping.
- **Step 5 partial circle slow**: Runs sequentially per sample. Large FASTQs (>10GB) take 10-30 min each.

## Batch Metadata

Batch TSVs are in `data/batches/`. 95 batches (batch_000 to batch_094), ~993 samples each, ~50,000 total samples.

To create new batch TSVs from an SRA search:
```bash
python scripts/search_sra.py --query "..." --output sra_results.tsv
python scripts/split_batches.py --input sra_results.tsv --output-dir data/batches/ --batch-size 993
```
