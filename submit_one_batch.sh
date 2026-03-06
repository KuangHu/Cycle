#!/bin/bash
# Process one batch end-to-end using N parallel SLURM nodes.
#
# Splits samples across N nodes for per-sample steps, consolidates for
# batch-level steps, then cleans up intermediate files.
#
# Usage:
#   bash submit_one_batch.sh BATCH_NUM [--nodes N]
#
# Example:
#   bash submit_one_batch.sh 6              # 5 nodes (default)
#   bash submit_one_batch.sh 6 --nodes 3    # 3 nodes
#
# Phases (dependency chain):
#
#   Phase 1 chunk_0 → Phase 2 chunk_0 ─┐
#   Phase 1 chunk_1 → Phase 2 chunk_1 ─┤
#   Phase 1 chunk_2 → Phase 2 chunk_2 ─┼→ Phase 3 → Phase 4 chunk_0 ─┐
#   Phase 1 chunk_3 → Phase 2 chunk_3 ─┤            Phase 4 chunk_1 ─┤
#   Phase 1 chunk_4 → Phase 2 chunk_4 ─┘            Phase 4 chunk_2 ─┼→ Phase 5
#                                                    Phase 4 chunk_3 ─┤
#                                                    Phase 4 chunk_4 ─┘
#
#   1  Download + Resolve + Align + Sniffles + Format   (N nodes, ~2-3 days)
#   2  ORF annotation + Guide finder                    (N nodes, ~1 hour)
#   3  Clustering + Novelty                             (1 node, ~6 hours)
#   4  Partial circle                                   (N nodes, ~2 hours)
#   5  Cleanup                                          (1 node, ~1 hour)
#
# All N nodes start immediately in Phase 1. Each chunk is independent.
# Resolve/index is idempotent — multiple nodes can do it safely.
#
# Automatically detects samples with existing guide.json and skips
# them for phases 1-2 (but still downloads their FASTQs for phase 4).
# Cleanup deletes FASTQs, BAMs, VCFs, circle intermediates, and assembly dirs.

set -e

# ── Parse arguments ────────────────────────────────────────────────
BATCH_NUM=""
NODES=5
while [[ $# -gt 0 ]]; do
    case $1 in
        --nodes) NODES=$2; shift 2 ;;
        *)
            if [[ -z "$BATCH_NUM" ]]; then
                BATCH_NUM=$1; shift
            else
                echo "Unknown argument: $1"; exit 1
            fi
            ;;
    esac
done

if [[ -z "$BATCH_NUM" ]]; then
    echo "Usage: bash submit_one_batch.sh BATCH_NUM [--nodes N]"
    echo ""
    echo "Phases:"
    echo "  1  Download + Resolve + Align + Format   (N nodes, ~2-3 days, start immediately)"
    echo "  2  ORF annotation + Guide finder         (N nodes, ~1 hour)"
    echo "  3  Clustering + Novelty                  (1 node, ~6 hours)"
    echo "  4  Partial circle                        (N nodes, ~2 hours)"
    echo "  5  Cleanup                               (1 node, ~1 hour)"
    exit 1
fi

batch_name=$(printf 'batch_%03d' "$BATCH_NUM")

# ── Paths ──────────────────────────────────────────────────────────
CYCLE_DIR=/home/kuangh/tools/Cycle
OUTROOT=/groups/rubin/projects/kuang/out/IS_cycle
BATCH_DIR=$OUTROOT/$batch_name
FMT_DIR=$BATCH_DIR/is_formatter_output
METADATA=$CYCLE_DIR/data/batches/${batch_name}.tsv
LOGDIR=$HOME/logs
CONDA_ENV=opfi
ENV_SETUP="export PATH=/home/kuangh/miniconda3/envs/${CONDA_ENV}/bin:\$PATH"

PIPELINE=$CYCLE_DIR/scripts/run_pipeline.py
ORF_ANNOTATOR=$CYCLE_DIR/scripts/run_orf_annotator.py
GUIDE_FINDER=$CYCLE_DIR/scripts/run_guide_finder.py
CLUSTERING=$CYCLE_DIR/scripts/run_system_clustering.py
NOVELTY=$CYCLE_DIR/scripts/run_novelty_annotator.py
PARTIAL_CIRCLE=$CYCLE_DIR/scripts/run_partial_circle.py
ISFINDER=$OUTROOT/batch_000/is_reference/ISfinder_raw.fna

mkdir -p "$LOGDIR"

if [[ ! -f "$METADATA" ]]; then
    echo "ERROR: Metadata not found: $METADATA"
    exit 1
fi

# ── Split metadata into chunks ─────────────────────────────────────
CHUNK_DIR=$BATCH_DIR/chunks
mkdir -p "$CHUNK_DIR"

echo "=========================================="
echo "Batch: $batch_name  Nodes: $NODES"
echo "=========================================="
echo ""
echo "Splitting metadata into $NODES chunks..."

python3 -c "
import pandas as pd
import os

metadata = pd.read_csv('$METADATA', sep='\t')
n_nodes = $NODES
fmt_dir = '$FMT_DIR'
chunk_dir = '$CHUNK_DIR'

# Detect completed samples (have guide.json)
completed = set()
if os.path.isdir(fmt_dir):
    for d in os.listdir(fmt_dir):
        guide = os.path.join(fmt_dir, d, f'{d}_is_records_guide.json')
        if os.path.exists(guide):
            completed.add(d)

print(f'Total samples: {len(metadata)}')
print(f'Already completed (have guide.json): {len(completed)}')
print(f'Incomplete: {len(metadata) - len(completed)}')
print()

# Split ALL samples into N chunks (same split for all phases)
chunk_size = (len(metadata) + n_nodes - 1) // n_nodes
for i in range(n_nodes):
    chunk = metadata.iloc[i*chunk_size:(i+1)*chunk_size]
    if len(chunk) == 0:
        continue

    # All samples in chunk (for download + partial circle)
    all_path = os.path.join(chunk_dir, f'chunk_{i}_all.tsv')
    chunk.to_csv(all_path, sep='\t', index=False)

    # Incomplete samples in chunk (for align + format + orf + guide)
    incomplete = chunk[~chunk['srr_accession'].isin(completed)]
    inc_path = os.path.join(chunk_dir, f'chunk_{i}_incomplete.tsv')
    incomplete.to_csv(inc_path, sep='\t', index=False)

    print(f'  chunk_{i}: {len(chunk)} total, {len(incomplete)} incomplete')
"

echo ""

# ── Phase 1: Download + Resolve + Align + Sniffles + Format (N nodes) ──
# All N nodes start immediately — no dependency bottleneck.
# Resolve/index is idempotent: first node to finish downloads refs, others skip.
echo "Phase 1: Download + Resolve + Align + Sniffles + Format ($NODES nodes)"

P1_JOBIDS=()
for i in $(seq 0 $((NODES - 1))); do
    ALL_TSV="$CHUNK_DIR/chunk_${i}_all.tsv"
    INC_TSV="$CHUNK_DIR/chunk_${i}_incomplete.tsv"
    CDIR="$CHUNK_DIR/chunk_${i}"

    [[ ! -f "$ALL_TSV" ]] && continue

    n_incomplete=$(tail -n +2 "$INC_TSV" 2>/dev/null | wc -l)

    # Download ALL samples in chunk (needed for partial circle later)
    # Resolve + Index references (idempotent across nodes)
    P1_CMD="$ENV_SETUP && \\
echo '=== Phase 1 chunk $i: Download + Resolve + Index ===' && \\
python $PIPELINE \\
  --metadata '$ALL_TSV' \\
  --fastq-dir '$CDIR/sra_downloads' \\
  --ref-dir '$BATCH_DIR/reference_genomes' \\
  --align-dir '$CDIR/alignments' \\
  --steps download resolve index \\
  --threads 44"

    # Align + Sniffles + Circle + Format for INCOMPLETE samples only
    if [[ "$n_incomplete" -gt 0 ]]; then
        P1_CMD="$P1_CMD && \\
echo '=== Phase 1 chunk $i: Align ($n_incomplete samples) ===' && \\
python $PIPELINE \\
  --metadata '$INC_TSV' \\
  --fastq-dir '$CDIR/sra_downloads' \\
  --ref-dir '$BATCH_DIR/reference_genomes' \\
  --align-dir '$CDIR/alignments' \\
  --steps align \\
  --threads 44 && \\
echo '=== Phase 1 chunk $i: Sniffles + Circle + Format ===' && \\
python $PIPELINE \\
  --metadata '$INC_TSV' \\
  --fastq-dir '$CDIR/sra_downloads' \\
  --ref-dir '$BATCH_DIR/reference_genomes' \\
  --align-dir '$CDIR/alignments' \\
  --sniffles-dir '$CDIR/sniffles_output' \\
  --circle-dir '$CDIR/circle_output' \\
  --formatter-dir '$FMT_DIR' \\
  --steps sniffles circle format \\
  --parallel 6 --threads 8"
    fi

    P1_CMD="$P1_CMD && echo '=== Phase 1 chunk $i done ==='"

    jobid=$(sbatch --parsable \
        --job-name="${batch_name}_p1_${i}" \
        --partition=standard --qos=standard \
        --nodes=1 --cpus-per-task=48 --mem=192G \
        --time=3-00:00:00 \
        --output="$LOGDIR/${batch_name}_p1_${i}_%j.out" \
        --error="$LOGDIR/${batch_name}_p1_${i}_%j.err" \
        --wrap="$P1_CMD")

    P1_JOBIDS+=("$jobid")
    echo "  chunk $i: Job $jobid (${n_incomplete} incomplete + download all)"
done

P1_DEP=$(IFS=:; echo "${P1_JOBIDS[*]}")

# ── Phase 2: ORF + Guide (N nodes, incomplete only) ────────────────
# Each chunk depends only on its own Phase 1 job.
echo ""
echo "Phase 2: ORF annotation + Guide finder ($NODES nodes)"

P2_JOBIDS=()
for i in $(seq 0 $((NODES - 1))); do
    INC_TSV="$CHUNK_DIR/chunk_${i}_incomplete.tsv"
    [[ ! -f "$INC_TSV" ]] && continue
    n_incomplete=$(tail -n +2 "$INC_TSV" 2>/dev/null | wc -l)
    [[ "$n_incomplete" -eq 0 ]] && continue

    P2_CMD="$ENV_SETUP && \\
echo '=== Phase 2 chunk $i: ORF annotation ===' && \\
python $ORF_ANNOTATOR \\
  --input-dir '$FMT_DIR' \\
  --sample-list '$INC_TSV' \\
  --parallel 20 && \\
echo '=== Phase 2 chunk $i: Guide finder ===' && \\
python $GUIDE_FINDER \\
  --input-dir '$FMT_DIR' \\
  --sample-list '$INC_TSV' \\
  --parallel 20 && \\
echo '=== Phase 2 chunk $i done ==='"

    jobid=$(sbatch --parsable \
        --job-name="${batch_name}_p2_${i}" \
        --partition=standard --qos=standard \
        --nodes=1 --cpus-per-task=48 --mem=192G \
        --time=6:00:00 \
        --output="$LOGDIR/${batch_name}_p2_${i}_%j.out" \
        --error="$LOGDIR/${batch_name}_p2_${i}_%j.err" \
        --dependency="afterok:${P1_JOBIDS[$i]}" \
        --wrap="$P2_CMD")

    P2_JOBIDS+=("$jobid")
    echo "  chunk $i: Job $jobid (dep: ${P1_JOBIDS[$i]})"
done

if [[ ${#P2_JOBIDS[@]} -eq 0 ]]; then
    echo "  (no incomplete samples — skipping phase 2)"
    P2_DEP="$P1_DEP"
else
    P2_DEP=$(IFS=:; echo "${P2_JOBIDS[*]}")
fi

# ── Phase 3: Clustering + Novelty (1 node) ─────────────────────────
# Depends on ALL Phase 2 jobs (needs all guide.json files).
echo ""
echo "Phase 3: Clustering + Novelty (1 node)"

CLUSTER_DIR=$OUTROOT/system_clustering_${batch_name}
NOVELTY_DIR=$OUTROOT/novelty_${batch_name}

# Determine ISfinder path
if [[ -f "$BATCH_DIR/is_reference/ISfinder_raw.fna" ]]; then
    ISFINDER_PATH="$BATCH_DIR/is_reference/ISfinder_raw.fna"
else
    ISFINDER_PATH="$ISFINDER"
fi

P3_CMD="$ENV_SETUP && \\
echo '=== Phase 3: System clustering ===' && \\
python $CLUSTERING \\
  --input-dirs '$FMT_DIR' \\
  --output-dir '$CLUSTER_DIR' \\
  --threads 44 && \\
echo '=== Phase 3: Novelty annotation ===' && \\
python $NOVELTY \\
  --input-dirs '$FMT_DIR' \\
  --clusters '$CLUSTER_DIR/system_clusters.json' \\
  --output-dir '$NOVELTY_DIR' \\
  --isfinder-fasta '$ISFINDER_PATH' \\
  --threads 44 && \\
echo '=== Phase 3 done ==='"

P3_JOBID=$(sbatch --parsable \
    --job-name="${batch_name}_p3" \
    --partition=standard --qos=standard \
    --nodes=1 --cpus-per-task=48 --mem=192G \
    --time=12:00:00 \
    --output="$LOGDIR/${batch_name}_p3_%j.out" \
    --error="$LOGDIR/${batch_name}_p3_%j.err" \
    --dependency="afterok:${P2_DEP}" \
    --wrap="$P3_CMD")

echo "  Job $P3_JOBID (dep: all Phase 2)"

# ── Phase 4: Partial circle (N nodes, all samples) ─────────────────
# Depends on Phase 3 (needs clustering done).
# Uses FASTQs from Phase 1 (still in chunks/chunk_N/sra_downloads/).
echo ""
echo "Phase 4: Partial circle ($NODES nodes)"

PC_DIR=$BATCH_DIR/partial_circle_output

P4_JOBIDS=()
for i in $(seq 0 $((NODES - 1))); do
    ALL_TSV="$CHUNK_DIR/chunk_${i}_all.tsv"
    CDIR="$CHUNK_DIR/chunk_${i}"
    [[ ! -f "$ALL_TSV" ]] && continue

    # Write per-chunk Python script (avoids nested quoting issues)
    PC_SCRIPT="$CHUNK_DIR/run_partial_circle_chunk_${i}.py"
    cat > "$PC_SCRIPT" << PCEOF
import os, subprocess, sys
sys.path.insert(0, "$CYCLE_DIR")
from cycle.utils import find_fastq
import pandas as pd

chunk_meta = pd.read_csv("$ALL_TSV", sep="\t")
sra_dir = "$CDIR/sra_downloads"
pc_dir = "$PC_DIR"
fmt_dir = "$FMT_DIR"
os.makedirs(pc_dir, exist_ok=True)

samples = []
for _, row in chunk_meta.iterrows():
    sid = row["srr_accession"]
    if os.path.isdir(os.path.join(pc_dir, sid)):
        continue
    jf = os.path.join(fmt_dir, sid, f"{sid}_is_records_guide.json")
    if not os.path.exists(jf):
        continue
    fq = find_fastq(sra_dir, sid)
    if fq:
        samples.append((sid, jf, str(fq)))

print(f"Partial circle chunk $i: {len(samples)} samples to process")
for sid, jf, fq in samples:
    out = os.path.join(pc_dir, sid)
    cmd = [
        "python", "$PARTIAL_CIRCLE",
        "--is-records", jf,
        "--fastq", fq,
        "--sample-id", sid,
        "--output-dir", out,
        "--threads", "8",
    ]
    print(f"  Running {sid}...")
    subprocess.run(cmd, check=False)
print(f"Partial circle chunk $i: done")
PCEOF

    P4_CMD="$ENV_SETUP && \\
echo '=== Phase 4 chunk $i: Partial circle ===' && \\
python3 '$PC_SCRIPT' && \\
echo '=== Phase 4 chunk $i done ==='"

    jobid=$(sbatch --parsable \
        --job-name="${batch_name}_p4_${i}" \
        --partition=standard --qos=standard \
        --nodes=1 --cpus-per-task=48 --mem=192G \
        --time=12:00:00 \
        --output="$LOGDIR/${batch_name}_p4_${i}_%j.out" \
        --error="$LOGDIR/${batch_name}_p4_${i}_%j.err" \
        --dependency="afterok:${P3_JOBID}" \
        --wrap="$P4_CMD")

    P4_JOBIDS+=("$jobid")
    echo "  chunk $i: Job $jobid"
done

P4_DEP=$(IFS=:; echo "${P4_JOBIDS[*]}")

# ── Phase 5: Cleanup (1 node) ──────────────────────────────────────
# Deletes chunk working dirs (FASTQs, BAMs, VCFs) and assembly dirs.
echo ""
echo "Phase 5: Cleanup (1 node)"

P5_CMD="$ENV_SETUP && \\
echo '=== Phase 5: Cleanup ===' && \\
echo 'Deleting chunk working directories (FASTQs, BAMs, VCFs, circle intermediates)...' && \\
rm -rf '$CHUNK_DIR' && \\
echo 'Deleting assembly directories...' && \\
find '$FMT_DIR' -maxdepth 2 -type d -name assembly -exec rm -rf {} + 2>/dev/null; \\
echo '=== Cleanup done ===' && \\
du -sh '$BATCH_DIR'"

P5_JOBID=$(sbatch --parsable \
    --job-name="${batch_name}_p5" \
    --partition=standard --qos=standard \
    --nodes=1 --cpus-per-task=4 --mem=16G \
    --time=6:00:00 \
    --output="$LOGDIR/${batch_name}_p5_%j.out" \
    --error="$LOGDIR/${batch_name}_p5_%j.err" \
    --dependency="afterok:${P4_DEP}" \
    --wrap="$P5_CMD")

echo "  Job $P5_JOBID"

# ── Summary ────────────────────────────────────────────────────────
TOTAL_JOBS=$((${#P1_JOBIDS[@]} + ${#P2_JOBIDS[@]} + 1 + ${#P4_JOBIDS[@]} + 1))
echo ""
echo "=========================================="
echo "Submitted $TOTAL_JOBS jobs for $batch_name"
echo "=========================================="
echo ""
echo "Phase 1 (download+format): ${P1_JOBIDS[*]}  ← start immediately"
echo "Phase 2 (orf+guide):       ${P2_JOBIDS[*]:-skipped}"
echo "Phase 3 (cluster+novelty): $P3_JOBID"
echo "Phase 4 (partial circle):  ${P4_JOBIDS[*]}"
echo "Phase 5 (cleanup):         $P5_JOBID"
echo ""
echo "Monitor: squeue -u \$(whoami) -n ${batch_name}_p"
echo "Logs:    tail -f $LOGDIR/${batch_name}_p*.err"
