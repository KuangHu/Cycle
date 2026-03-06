#!/bin/bash
# Submit the full batch pipeline in 5 chained steps.
#
# Usage:
#   bash submit_all.sh --step 1 [start] [end]    # Download + Resolve + Index
#   bash submit_all.sh --step 2 [start] [end]    # Align
#   bash submit_all.sh --step 3 [start] [end]    # Sniffles + Circle + Format
#   bash submit_all.sh --step 4 [start] [end]    # ORF annotation + Guide finder
#   bash submit_all.sh --step 5 [start] [end]    # System clustering + Novelty + Partial circle
#
# Chain steps with --dep:
#   bash submit_all.sh --step 2 --dep 1 21 25
#   bash submit_all.sh --step 3 --dep 2 21 25
#   bash submit_all.sh --step 4 --dep 3 21 25
#   bash submit_all.sh --step 5 --dep 4 21 25
#
# Full end-to-end (5 batches):
#   for s in 1 2 3 4 5; do
#     dep=""; [ $s -gt 1 ] && dep="--dep $((s-1))"
#     bash submit_all.sh --step $s $dep 21 25
#   done
#
# Steps overview:
#   1  Download + Resolve + Index       (whole node, ~2 days)
#   2  Align                            (whole node, ~2 days)
#   3  Sniffles + Circle + Format       (whole node, ~2 days)
#   4  ORF annotation + Guide finder    (whole node, ~6 hours)
#   5  Clustering + Novelty + Partial   (whole node, ~12 hours)

set -e

# ── Parse arguments ──────────────────────────────────────────────────
STEP=""
DEP_STEP=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --step) STEP="$2"; shift 2 ;;
        --dep)  DEP_STEP="$2"; shift 2 ;;
        *) break ;;
    esac
done

if [ -z "$STEP" ] || ! [[ "$STEP" =~ ^[12345]$ ]]; then
    echo "Usage: bash submit_all.sh --step {1|2|3|4|5} [--dep {1|2|3|4|5}] [start] [end]"
    echo ""
    echo "  Step 1: Download + Resolve + Index       (48 CPUs, 192G, 2 days)"
    echo "  Step 2: Align                            (48 CPUs, 192G, 2 days)"
    echo "  Step 3: Sniffles + Circle + Format       (48 CPUs, 192G, 2 days)"
    echo "  Step 4: ORF annotation + Guide finder    (48 CPUs, 192G, 6 hours)"
    echo "  Step 5: Clustering + Novelty + Partial   (48 CPUs, 192G, 12 hours)"
    echo ""
    echo "  --dep N: wait for step N to finish (SLURM afterok dependency)"
    exit 1
fi

START=${1:-5}
END=${2:-94}

# ── Paths ────────────────────────────────────────────────────────────
CYCLE_DIR=/home/kuangh/tools/Cycle
PIPELINE=$CYCLE_DIR/scripts/run_pipeline.py
PREPARE_META=$CYCLE_DIR/scripts/prepare_batch_metadata.py
GUIDE_FINDER=$CYCLE_DIR/scripts/run_guide_finder.py
ORF_ANNOTATOR=$CYCLE_DIR/scripts/run_orf_annotator.py
CLUSTERING=$CYCLE_DIR/scripts/run_system_clustering.py
NOVELTY=$CYCLE_DIR/scripts/run_novelty_annotator.py
PARTIAL_CIRCLE=$CYCLE_DIR/scripts/run_partial_circle.py

OUTROOT=/groups/rubin/projects/kuang/out/IS_cycle
LOGDIR=$HOME/logs
BATCHDIR=$CYCLE_DIR/data/batches
ISFINDER=$OUTROOT/batch_000/is_reference/ISfinder_raw.fna
CONDA_ENV=opfi

mkdir -p "$LOGDIR"

# ── Per-step resources ───────────────────────────────────────────────
case $STEP in
    1) CPUS=48; MEM="192G"; TIME="2-00:00:00" ;;
    2) CPUS=48; MEM="192G"; TIME="2-00:00:00" ;;
    3) CPUS=48; MEM="192G"; TIME="2-00:00:00" ;;
    4) CPUS=48; MEM="192G"; TIME="6:00:00"    ;;
    5) CPUS=48; MEM="192G"; TIME="12:00:00"   ;;
esac

echo "=========================================="
echo "Step: $STEP  (CPUs=$CPUS, Mem=$MEM, Time=$TIME)"
[ -n "$DEP_STEP" ] && echo "Dependency: afterok step $DEP_STEP"
echo "Batches: $(printf '%03d' $START) to $(printf '%03d' $END)"
echo "=========================================="
echo ""

count=0
for i in $(seq $START $END); do
    batch_name=$(printf 'batch_%03d' $i)
    metadata="$BATCHDIR/${batch_name}.tsv"
    batch_dir="$OUTROOT/$batch_name"
    fmt_dir="$batch_dir/is_formatter_output"

    # ── Resolve dependency ───────────────────────────────────────────
    DEP_FLAG=""
    if [ -n "$DEP_STEP" ]; then
        dep_jobid=$(squeue -u "$(whoami)" -n "${batch_name}_s${DEP_STEP}" -h -o "%i" 2>/dev/null | head -1)
        if [ -n "$dep_jobid" ]; then
            DEP_FLAG="--dependency=afterok:${dep_jobid}"
            echo "  $batch_name: depends on job $dep_jobid (step $DEP_STEP)"
        else
            echo "  $batch_name: no running step $DEP_STEP job found, submitting immediately"
        fi
    fi

    # ── Pre-flight checks (skip if dependency set) ───────────────────
    if [ -z "$DEP_FLAG" ]; then
        case $STEP in
            1) [ ! -f "$metadata" ] && echo "  SKIP $batch_name — metadata not found" && continue ;;
            2) [ ! -d "$batch_dir/sra_downloads" ] && echo "  SKIP $batch_name — no sra_downloads/" && continue ;;
            3) [ ! -d "$batch_dir/alignments" ] && echo "  SKIP $batch_name — no alignments/" && continue ;;
            4) [ ! -d "$fmt_dir" ] && echo "  SKIP $batch_name — no is_formatter_output/" && continue ;;
            5) [ ! -d "$fmt_dir" ] && echo "  SKIP $batch_name — no is_formatter_output/" && continue ;;
        esac
    fi

    # ── Build command ────────────────────────────────────────────────
    CMD="export PATH=/home/kuangh/miniconda3/envs/${CONDA_ENV}/bin:\$PATH"

    case $STEP in
        1)
            CMD="$CMD && \\
echo '=== Step 1: Download + Resolve + Index ===' && \\
python $PIPELINE \\
  --metadata '$metadata' \\
  --outdir '$batch_dir' \\
  --steps download resolve index \\
  --threads 44"
            ;;
        2)
            CMD="$CMD && \\
echo '=== Step 2: Align ===' && \\
python $PIPELINE \\
  --metadata '$metadata' \\
  --outdir '$batch_dir' \\
  --steps align \\
  --threads 44 \\
  --sort-memory 4G"
            ;;
        3)
            CMD="$CMD && \\
echo '=== Step 3: Sniffles + Circle + Format ===' && \\
python $PREPARE_META '$batch_dir' && \\
python $PIPELINE \\
  --metadata '$batch_dir/metadata_for_sniffles.tsv' \\
  --outdir '$batch_dir' \\
  --steps sniffles circle format \\
  --parallel 6 \\
  --threads 8"
            ;;
        4)
            CMD="$CMD && \\
echo '=== Step 4a: ORF annotation ===' && \\
python $ORF_ANNOTATOR \\
  --input-dir '$fmt_dir' \\
  --parallel 20 && \\
echo '=== Step 4b: Guide finder ===' && \\
python $GUIDE_FINDER \\
  --input-dir '$fmt_dir' \\
  --parallel 20"
            ;;
        5)
            # Clustering + Novelty
            cluster_dir="$OUTROOT/system_clustering_${batch_name}"
            novelty_dir="$OUTROOT/novelty_${batch_name}"
            pc_dir="$batch_dir/partial_circle_output"

            # Determine ISfinder path
            if [ -f "$batch_dir/is_reference/ISfinder_raw.fna" ]; then
                isfinder="$batch_dir/is_reference/ISfinder_raw.fna"
            else
                isfinder="$ISFINDER"
            fi

            CMD="$CMD && \\
echo '=== Step 5a: System clustering ===' && \\
python $CLUSTERING \\
  --input-dirs '$fmt_dir' \\
  --output-dir '$cluster_dir' \\
  --threads 44 && \\
echo '=== Step 5b: Novelty annotation ===' && \\
python $NOVELTY \\
  --input-dirs '$fmt_dir' \\
  --clusters '$cluster_dir/system_clusters.json' \\
  --output-dir '$novelty_dir' \\
  --isfinder-fasta '$isfinder' \\
  --threads 44 && \\
echo '=== Step 5c: Partial circle detection ===' && \\
python -c \"
import json, os, glob, subprocess, sys
sys.path.insert(0, '$CYCLE_DIR')
from cycle.utils import find_fastq

fmt_dir = '$fmt_dir'
sra_dir = '$batch_dir/sra_downloads'
pc_dir = '$pc_dir'
os.makedirs(pc_dir, exist_ok=True)

# Find all samples with guide JSONs
samples = []
for jf in sorted(glob.glob(os.path.join(fmt_dir, '*', '*_is_records_guide.json'))):
    sid = os.path.basename(os.path.dirname(jf))
    if os.path.isdir(os.path.join(pc_dir, sid)):
        continue  # already done
    fq = find_fastq(sra_dir, sid)
    if fq:
        samples.append((sid, jf, str(fq)))

print(f'Partial circle: {len(samples)} samples to process')
for sid, jf, fq in samples:
    out = os.path.join(pc_dir, sid)
    cmd = [
        'python', '$PARTIAL_CIRCLE',
        '--is-records', jf,
        '--fastq', fq,
        '--sample-id', sid,
        '--output-dir', out,
        '--threads', '8',
    ]
    print(f'  Running {sid}...')
    subprocess.run(cmd, check=False)
print('Partial circle: done')
\""
            ;;
    esac

    CMD="$CMD && echo '=== Done ==='"

    # ── Submit ───────────────────────────────────────────────────────
    jobid=$(sbatch --parsable \
        --job-name="${batch_name}_s${STEP}" \
        --partition=standard \
        --qos=standard \
        --nodes=1 \
        --cpus-per-task=$CPUS \
        --mem=$MEM \
        --time=$TIME \
        --output="$LOGDIR/${batch_name}_step${STEP}_%j.out" \
        --error="$LOGDIR/${batch_name}_step${STEP}_%j.err" \
        $DEP_FLAG \
        --wrap="$CMD")

    echo "  $batch_name: Job $jobid"
    count=$((count + 1))
done

echo ""
echo "=========================================="
echo "Submitted $count jobs (step=$STEP, cpus=$CPUS, mem=$MEM)"
echo "=========================================="
echo ""
echo "Monitor: squeue -u \$(whoami)"
echo "Logs:    tail -f $LOGDIR/batch_*_step${STEP}_*.err"
