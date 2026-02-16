#!/bin/bash
# Run Sniffles2 + circle detection on batches 0-4
# Each batch runs independently

# Create logs directory
mkdir -p /home/kuangh/logs

# Batch 0
echo "Submitting batch_000..."
sbatch --job-name=snif_circ_000 \
       --partition=standard \
       --qos=standard \
       --nodes=1 \
       --cpus-per-task=48 \
       --mem=180G \
       --time=2-00:00:00 \
       --output=/home/kuangh/logs/batch_000_sniffles_circle_%j.out \
       --error=/home/kuangh/logs/batch_000_sniffles_circle_%j.err \
       --wrap="
eval \"\$(conda shell.bash hook)\" && conda activate opfi && \
cd /groups/rubin/projects/kuang/out/IS_cycle/batch_000 && \
echo 'Starting batch_000 at \$(date)' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --outdir /groups/rubin/projects/kuang/out/IS_cycle/batch_000 \
  --steps sniffles circle \
  --sniffles-parallel 6 \
  --threads 8 && \
echo 'Completed batch_000 at \$(date)'
"

# Batch 1
echo "Submitting batch_001..."
sbatch --job-name=snif_circ_001 \
       --partition=standard \
       --qos=standard \
       --nodes=1 \
       --cpus-per-task=48 \
       --mem=180G \
       --time=2-00:00:00 \
       --output=/home/kuangh/logs/batch_001_sniffles_circle_%j.out \
       --error=/home/kuangh/logs/batch_001_sniffles_circle_%j.err \
       --wrap="
eval \"\$(conda shell.bash hook)\" && conda activate opfi && \
cd /groups/rubin/projects/kuang/out/IS_cycle/batch_001 && \
echo 'Starting batch_001 at \$(date)' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --outdir /groups/rubin/projects/kuang/out/IS_cycle/batch_001 \
  --steps sniffles circle \
  --sniffles-parallel 6 \
  --threads 8 && \
echo 'Completed batch_001 at \$(date)'
"

# Batch 2
echo "Submitting batch_002..."
sbatch --job-name=snif_circ_002 \
       --partition=standard \
       --qos=standard \
       --nodes=1 \
       --cpus-per-task=48 \
       --mem=180G \
       --time=2-00:00:00 \
       --output=/home/kuangh/logs/batch_002_sniffles_circle_%j.out \
       --error=/home/kuangh/logs/batch_002_sniffles_circle_%j.err \
       --wrap="
eval \"\$(conda shell.bash hook)\" && conda activate opfi && \
cd /groups/rubin/projects/kuang/out/IS_cycle/batch_002 && \
echo 'Starting batch_002 at \$(date)' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --outdir /groups/rubin/projects/kuang/out/IS_cycle/batch_002 \
  --steps sniffles circle \
  --sniffles-parallel 6 \
  --threads 8 && \
echo 'Completed batch_002 at \$(date)'
"

# Batch 3
echo "Submitting batch_003..."
sbatch --job-name=snif_circ_003 \
       --partition=standard \
       --qos=standard \
       --nodes=1 \
       --cpus-per-task=48 \
       --mem=180G \
       --time=2-00:00:00 \
       --output=/home/kuangh/logs/batch_003_sniffles_circle_%j.out \
       --error=/home/kuangh/logs/batch_003_sniffles_circle_%j.err \
       --wrap="
eval \"\$(conda shell.bash hook)\" && conda activate opfi && \
cd /groups/rubin/projects/kuang/out/IS_cycle/batch_003 && \
echo 'Starting batch_003 at \$(date)' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --outdir /groups/rubin/projects/kuang/out/IS_cycle/batch_003 \
  --steps sniffles circle \
  --sniffles-parallel 6 \
  --threads 8 && \
echo 'Completed batch_003 at \$(date)'
"

# Batch 4
echo "Submitting batch_004..."
sbatch --job-name=snif_circ_004 \
       --partition=standard \
       --qos=standard \
       --nodes=1 \
       --cpus-per-task=48 \
       --mem=180G \
       --time=2-00:00:00 \
       --output=/home/kuangh/logs/batch_004_sniffles_circle_%j.out \
       --error=/home/kuangh/logs/batch_004_sniffles_circle_%j.err \
       --wrap="
eval \"\$(conda shell.bash hook)\" && conda activate opfi && \
cd /groups/rubin/projects/kuang/out/IS_cycle/batch_004 && \
echo 'Starting batch_004 at \$(date)' && \
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
  --outdir /groups/rubin/projects/kuang/out/IS_cycle/batch_004 \
  --steps sniffles circle \
  --sniffles-parallel 6 \
  --threads 8 && \
echo 'Completed batch_004 at \$(date)'
"

echo ""
echo "All 5 batches submitted!"
echo "Check job status: squeue -u $USER"
echo "Monitor logs: tail -f /home/kuangh/logs/batch_00*_sniffles_circle_*.out"
