#!/bin/bash
# Submit Sniffles2 + circle detection for batches 0-4
# Run this script: bash submit_batches_0to4.sh

set -e

# Create logs directory
mkdir -p /home/kuangh/logs

echo "=========================================="
echo "Preparing and submitting batches 0-4"
echo "=========================================="

# Generate metadata for all batches first
for i in {0..4}; do
    batch_dir="/groups/rubin/projects/kuang/out/IS_cycle/batch_00$i"
    echo "Generating metadata for batch_00$i..."
    python /home/kuangh/tools/Cycle/scripts/prepare_batch_metadata.py "$batch_dir"
done

echo ""
echo "Submitting SLURM jobs..."
echo ""

# Submit each batch as a separate job
for i in {0..4}; do
    batch_name="batch_00$i"
    batch_dir="/groups/rubin/projects/kuang/out/IS_cycle/$batch_name"
    metadata="$batch_dir/metadata_for_sniffles.tsv"

    jobid=$(sbatch --parsable \
        --job-name="snif_$batch_name" \
        --partition=standard \
        --qos=standard \
        --nodes=1 \
        --cpus-per-task=48 \
        --mem=180G \
        --time=2-00:00:00 \
        --output="/home/kuangh/logs/${batch_name}_sniffles_circle_%j.out" \
        --error="/home/kuangh/logs/${batch_name}_sniffles_circle_%j.err" \
        --wrap="
eval \"\$(conda shell.bash hook)\" && \
conda activate opfi && \
echo '========================================' && \
echo 'Starting $batch_name at \$(date)' && \
echo 'Metadata: $metadata' && \
echo 'Output: $batch_dir' && \
echo '========================================' && \
echo '' && \
echo 'Step 1: Sniffles2 IS detection (6 organisms in parallel)' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --metadata '$metadata' \
  --outdir '$batch_dir' \
  --steps sniffles \
  --sniffles-parallel 6 \
  --threads 8 && \
echo '' && \
echo 'Step 2: Circle detection' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --metadata '$metadata' \
  --outdir '$batch_dir' \
  --steps circle \
  --threads 8 && \
echo '' && \
echo '========================================' && \
echo 'Completed $batch_name at \$(date)' && \
echo '========================================'
")

    echo "  $batch_name: Job $jobid"
done

echo ""
echo "=========================================="
echo "All 5 batches submitted!"
echo "=========================================="
echo ""
echo "Monitor jobs:"
echo "  squeue -u $USER"
echo ""
echo "Watch logs:"
echo "  tail -f /home/kuangh/logs/batch_00*_sniffles_circle_*.out"
echo ""
echo "Check status:"
echo "  ls -lh /groups/rubin/projects/kuang/out/IS_cycle/batch_00*/sniffles_output/"
echo "  ls -lh /groups/rubin/projects/kuang/out/IS_cycle/batch_00*/circle_output/"
echo ""
