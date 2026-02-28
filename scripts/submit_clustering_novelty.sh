#!/bin/bash
# Submit system clustering + novelty annotation for a range of batches.
#
# Output follows the canonical data structure (see DATA_STRUCTURE.md):
#   - system_clustering_batch_NNN/  at IS_cycle root
#   - novelty_batch_NNN/            at IS_cycle root
#
# Usage:
#   bash scripts/submit_clustering_novelty.sh 0 4      # batches 000-004
#   bash scripts/submit_clustering_novelty.sh 6 10     # batches 006-010
#   bash scripts/submit_clustering_novelty.sh 11 15    # batches 011-015

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 START_BATCH END_BATCH"
    echo "  e.g.: $0 0 4   (batches 000-004)"
    exit 1
fi

START=$1
END=$2

CYCLE_DIR="/home/kuangh/tools/Cycle"
BASE="/groups/rubin/projects/kuang/out/IS_cycle"
LOGDIR="/home/kuangh/logs"
ISFINDER="$BASE/batch_000/is_reference/ISfinder_raw.fna"

mkdir -p "$LOGDIR"

echo "=========================================="
echo "System clustering + novelty annotation"
echo "Batches: $(printf '%03d' $START) to $(printf '%03d' $END)"
echo "=========================================="

for i in $(seq $START $END); do
    batch=$(printf "batch_%03d" $i)
    batch_dir="$BASE/$batch"
    fmt_dir="$batch_dir/is_formatter_output"
    cluster_dir="$BASE/system_clustering_${batch}"
    novelty_dir="$BASE/novelty_${batch}"

    # Validate formatter output exists
    if [ ! -d "$fmt_dir" ]; then
        echo "  SKIP $batch: $fmt_dir not found"
        continue
    fi

    # Decide what to submit
    need_clust=true
    need_novel=true

    if [ -f "$cluster_dir/system_clusters.json" ]; then
        need_clust=false
    fi
    if [ -f "$novelty_dir/cluster_novelty_summary.tsv" ]; then
        need_novel=false
    fi

    if [ "$need_clust" = false ] && [ "$need_novel" = false ]; then
        echo "  SKIP $batch: clustering + novelty already done"
        continue
    fi

    # Determine ISfinder FASTA path
    isfinder_flag=""
    if [ -f "$batch_dir/is_reference/ISfinder_raw.fna" ]; then
        isfinder_flag="--isfinder-fasta $batch_dir/is_reference/ISfinder_raw.fna"
    else
        isfinder_flag="--isfinder-fasta $ISFINDER"
    fi

    if [ "$need_clust" = true ]; then
        # Submit clustering
        clust_jobid=$(sbatch --parsable \
            --job-name="clust_${batch}" \
            --partition=standard \
            --qos=standard \
            --nodes=1 \
            --cpus-per-task=48 \
            --mem=180G \
            --time=6:00:00 \
            --output="$LOGDIR/${batch}_clustering_%j.out" \
            --error="$LOGDIR/${batch}_clustering_%j.err" \
            --wrap="
eval \"\$(conda shell.bash hook)\" && \
conda activate opfi && \
echo '=== System clustering: $batch ===' && \
echo 'Start: '\$(date) && \
python $CYCLE_DIR/scripts/run_system_clustering.py \
    --input-dirs $fmt_dir \
    --output-dir $cluster_dir \
    --threads 48 && \
echo 'Done: '\$(date)
")

        if [ "$need_novel" = true ]; then
            # Chain novelty after clustering
            nov_jobid=$(sbatch --parsable \
                --job-name="novel_${batch}" \
                --partition=standard \
                --qos=standard \
                --nodes=1 \
                --cpus-per-task=48 \
                --mem=180G \
                --time=12:00:00 \
                --dependency=afterok:${clust_jobid} \
                --output="$LOGDIR/${batch}_novelty_%j.out" \
                --error="$LOGDIR/${batch}_novelty_%j.err" \
                --wrap="
eval \"\$(conda shell.bash hook)\" && \
conda activate opfi && \
echo '=== Novelty annotation: $batch ===' && \
echo 'Start: '\$(date) && \
python $CYCLE_DIR/scripts/run_novelty_annotator.py \
    --input-dirs $fmt_dir \
    --clusters $cluster_dir/system_clusters.json \
    --output-dir $novelty_dir \
    $isfinder_flag \
    --threads 48 && \
echo 'Done: '\$(date)
")
            echo "  $batch: clustering=$clust_jobid -> novelty=$nov_jobid"
        else
            echo "  $batch: clustering=$clust_jobid (novelty already done)"
        fi
    else
        # Clustering done, only submit novelty
        nov_jobid=$(sbatch --parsable \
            --job-name="novel_${batch}" \
            --partition=standard \
            --qos=standard \
            --nodes=1 \
            --cpus-per-task=48 \
            --mem=180G \
            --time=12:00:00 \
            --output="$LOGDIR/${batch}_novelty_%j.out" \
            --error="$LOGDIR/${batch}_novelty_%j.err" \
            --wrap="
eval \"\$(conda shell.bash hook)\" && \
conda activate opfi && \
echo '=== Novelty annotation: $batch ===' && \
echo 'Start: '\$(date) && \
python $CYCLE_DIR/scripts/run_novelty_annotator.py \
    --input-dirs $fmt_dir \
    --clusters $cluster_dir/system_clusters.json \
    --output-dir $novelty_dir \
    $isfinder_flag \
    --threads 48 && \
echo 'Done: '\$(date)
")
        echo "  $batch: novelty=$nov_jobid (clustering already done)"
    fi
done

echo ""
echo "=========================================="
echo "Done. Monitor: squeue -u kuangh"
echo "=========================================="
