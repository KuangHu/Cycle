#!/bin/bash
# Submit ORF annotator + guide finder for batches.
#
# Usage:
#   bash submit_orf_guide.sh 1 4      # submit batch_001 through batch_004

set -e

START=${1:?Usage: bash submit_orf_guide.sh START END}
END=${2:?Usage: bash submit_orf_guide.sh START END}

OUTROOT=/groups/rubin/projects/kuang/out/IS_cycle
LOGDIR=$HOME/logs
CODEDIR=/home/kuangh/tools/Cycle

mkdir -p "$LOGDIR"

echo "=========================================="
echo "ORF annotator + Guide finder"
echo "Batches: $(printf '%03d' $START) to $(printf '%03d' $END)"
echo "=========================================="

count=0
for i in $(seq $START $END); do
    batch_name=$(printf 'batch_%03d' $i)
    formatter_dir="$OUTROOT/$batch_name/is_formatter_output"

    if [ ! -d "$formatter_dir" ]; then
        echo "  SKIP $batch_name — no is_formatter_output/"
        continue
    fi

    CMD="export PATH=/home/kuangh/miniconda3/envs/opfi/bin:\$PATH && \
echo '=== ORF annotator: $batch_name ===' && \
python $CODEDIR/scripts/run_orf_annotator.py \
  --input-dir '$formatter_dir' \
  --parallel 44 && \
echo '=== Guide finder: $batch_name ===' && \
python $CODEDIR/scripts/run_guide_finder.py \
  --input-dir '$formatter_dir' \
  --parallel 44 && \
echo '=== Done: $batch_name ==='"

    jobid=$(sbatch --parsable \
        --job-name="${batch_name}_orf_guide" \
        --partition=standard \
        --qos=standard \
        --nodes=1 \
        --cpus-per-task=48 \
        --mem=192G \
        --time=2-00:00:00 \
        --output="$LOGDIR/${batch_name}_orf_guide_%j.out" \
        --error="$LOGDIR/${batch_name}_orf_guide_%j.err" \
        --wrap="$CMD")

    echo "  $batch_name: Job $jobid ($formatter_dir)"
    count=$((count + 1))
done

echo ""
echo "=========================================="
echo "Submitted $count jobs (48 CPUs, 192G each)"
echo "=========================================="
echo "Monitor: squeue -u \$(whoami)"
echo "Logs:    tail -f $LOGDIR/batch_*_orf_guide_*.err"
