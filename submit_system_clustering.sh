#!/bin/bash
#SBATCH --job-name=sys_cluster
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/groups/rubin/projects/kuang/out/IS_cycle/system_clustering_batch000/slurm_%j.log

set -euo pipefail

export PATH=/home/kuangh/miniconda3/envs/opfi/bin:$PATH

python /home/kuangh/tools/Cycle/scripts/run_system_clustering.py \
    --input-dirs /groups/rubin/projects/kuang/out/IS_cycle/batch_000/is_formatter_output \
    --output-dir /groups/rubin/projects/kuang/out/IS_cycle/system_clustering_batch000 \
    --threads 8 \
    --no-skip
