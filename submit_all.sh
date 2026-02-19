#!/bin/bash
# Submit pipeline for batches in three steps with right-sized resources.
#
# Usage:
#   bash submit_all.sh --step 1 [start] [end]    # Download + Resolve + Index
#   bash submit_all.sh --step 2 [start] [end]    # Align
#   bash submit_all.sh --step 3 [start] [end]    # Sniffles + Circle
#
# Note: partition is EXCLUSIVE (whole node per job), so all steps request full node.
#
# Examples:
#   bash submit_all.sh --step 1 5 94             # Download+resolve+index batches 005-094
#   bash submit_all.sh --step 2 --dep 1 5 94     # Align after step 1 finishes
#   bash submit_all.sh --step 3 --dep 2 5 94     # Sniffles+circle after step 2 finishes

set -e

# Parse arguments
STEP=""
DEP_STEP=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --step) STEP="$2"; shift 2 ;;
        --dep)  DEP_STEP="$2"; shift 2 ;;
        *) break ;;
    esac
done

if [ -z "$STEP" ] || ! [[ "$STEP" =~ ^[123]$ ]]; then
    echo "Usage: bash submit_all.sh --step {1|2|3} [--dep {1|2|3}] [start] [end]"
    echo ""
    echo "  Step 1: Download + Resolve + Index  (48 CPUs, 192G)"
    echo "  Step 2: Align                       (48 CPUs, 192G)"
    echo "  Step 3: Sniffles + Circle           (48 CPUs, 192G)"
    echo ""
    echo "  --dep N: wait for step N to finish (uses SLURM afterok dependency)"
    exit 1
fi

START=${1:-5}
END=${2:-94}

PIPELINE=/home/kuangh/tools/Cycle/scripts/run_pipeline.py
PREPARE_META=/home/kuangh/tools/Cycle/scripts/prepare_batch_metadata.py
OUTROOT=/groups/rubin/projects/kuang/out/IS_cycle
LOGDIR=$HOME/logs
BATCHDIR=/home/kuangh/tools/Cycle/data/batches

mkdir -p "$LOGDIR"

# Per-step SLURM resources
case $STEP in
    1) CPUS=48; MEM="192G"; TIME="2-00:00:00" ;;
    2) CPUS=48; MEM="192G"; TIME="2-00:00:00" ;;
    3) CPUS=48; MEM="192G"; TIME="2-00:00:00" ;;
esac

echo "=========================================="
echo "Step: $STEP  (CPUs=$CPUS, Mem=$MEM)"
[ -n "$DEP_STEP" ] && echo "Dependency: afterok step $DEP_STEP"
echo "Batches: $(printf '%03d' $START) to $(printf '%03d' $END)"
echo "=========================================="
echo ""

count=0
for i in $(seq $START $END); do
    batch_name=$(printf 'batch_%03d' $i)
    metadata="$BATCHDIR/${batch_name}.tsv"
    batch_dir="$OUTROOT/$batch_name"

    # Resolve dependency job ID for this batch
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

    # Step 1: need metadata file
    if [ "$STEP" = "1" ]; then
        if [ ! -f "$metadata" ]; then
            echo "  SKIP $batch_name — metadata not found"
            continue
        fi
    fi

    # Step 2: need downloaded FASTQs (skip check if dependency set)
    if [ "$STEP" = "2" ] && [ -z "$DEP_FLAG" ]; then
        if [ ! -d "$batch_dir/sra_downloads" ]; then
            echo "  SKIP $batch_name — no sra_downloads/ (run step 1 first)"
            continue
        fi
    fi

    # Step 3: need alignments (skip check if dependency set)
    if [ "$STEP" = "3" ] && [ -z "$DEP_FLAG" ]; then
        if [ ! -d "$batch_dir/alignments" ]; then
            echo "  SKIP $batch_name — no alignments/ (run step 2 first)"
            continue
        fi
    fi

    # Build the command based on step
    CMD="export PATH=/home/kuangh/miniconda3/envs/opfi/bin:\$PATH"

    if [ "$STEP" = "1" ]; then
        CMD="$CMD && \
echo '=== Step 1: Download + Resolve + Index ===' && \
python $PIPELINE \
  --metadata '$metadata' \
  --outdir '$batch_dir' \
  --steps download resolve index \
  --threads 44"
    fi

    if [ "$STEP" = "2" ]; then
        CMD="$CMD && \
echo '=== Step 2: Align ===' && \
python $PIPELINE \
  --metadata '$metadata' \
  --outdir '$batch_dir' \
  --steps align \
  --threads 44 \
  --sort-memory 4G"
    fi

    if [ "$STEP" = "3" ]; then
        CMD="$CMD && \
echo '=== Step 3: Sniffles + Circle ===' && \
python $PREPARE_META '$batch_dir' && \
python $PIPELINE \
  --metadata '$batch_dir/metadata_for_sniffles.tsv' \
  --outdir '$batch_dir' \
  --steps sniffles \
  --sniffles-parallel 6 \
  --threads 8 && \
python $PIPELINE \
  --metadata '$batch_dir/metadata_for_sniffles.tsv' \
  --outdir '$batch_dir' \
  --steps circle \
  --circle-parallel 6 \
  --threads 8"
    fi

    CMD="$CMD && echo '=== Done ==='"

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
