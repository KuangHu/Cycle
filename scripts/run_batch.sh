#!/bin/bash
# Launcher script for SLURM: activates conda and runs the pipeline.
# Usage: run_batch.sh <batch_tsv> <outdir> [threads] [sort_memory]

set -euo pipefail

BATCH_TSV="$1"
OUTDIR="$2"
THREADS="${3:-44}"
SORT_MEMORY="${4:-4G}"

eval "$(conda shell.bash hook)"
conda activate opfi

python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
    --metadata "$BATCH_TSV" \
    --outdir "$OUTDIR" \
    --threads "$THREADS" \
    --sort-memory "$SORT_MEMORY"
