#!/bin/bash
#SBATCH --job-name=sniffles_circle
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/kuangh/logs/sniffles_circle_%A_%a.out
#SBATCH --error=/home/kuangh/logs/sniffles_circle_%A_%a.err

# Run Sniffles2 + circle detection on a single batch
# Usage: sbatch --array=0-4 run_sniffles_circle_batch.sh

set -e

# Batch directories
BATCHES=(
    /groups/rubin/projects/kuang/out/IS_cycle/batch_000
    /groups/rubin/projects/kuang/out/IS_cycle/batch_001
    /groups/rubin/projects/kuang/out/IS_cycle/batch_002
    /groups/rubin/projects/kuang/out/IS_cycle/batch_003
    /groups/rubin/projects/kuang/out/IS_cycle/batch_004
)

BATCH_DIR="${BATCHES[$SLURM_ARRAY_TASK_ID]}"
BATCH_NAME=$(basename "$BATCH_DIR")

echo "=========================================="
echo "Processing: $BATCH_NAME"
echo "Batch directory: $BATCH_DIR"
echo "Time: $(date)"
echo "=========================================="

# Create temporary metadata from alignments
TEMP_METADATA="/tmp/${BATCH_NAME}_metadata_${SLURM_JOB_ID}.tsv"

cd "$BATCH_DIR"

# Extract organism names from alignment directory structure
# Assuming alignments are named like: SRR*.sorted.bam
echo "Discovering samples from alignments..."
python3 << 'PYTHON'
import pandas as pd
from pathlib import Path
import sys
import subprocess

batch_dir = Path(sys.argv[1])
alignments = list((batch_dir / "alignments").glob("*.sorted.bam"))

# Get organism info from reference genomes
ref_dir = batch_dir / "reference_genomes"
organisms = {}
for ref_path in ref_dir.glob("*/"):
    accession = ref_path.name
    # Try to find organism name from fasta file
    fasta = list(ref_path.glob("*.fna"))
    if fasta:
        organisms[accession] = ref_path.name

print(f"Found {len(alignments)} alignments", file=sys.stderr)

# Create metadata DataFrame
records = []
for bam in alignments:
    # Extract accession from filename (e.g., SRR21465445.sorted.bam -> SRR21465445)
    srr = bam.stem.replace('.sorted', '')
    records.append({
        'srr_accession': srr,
        'organism': 'unknown',  # Will be determined from reference mapping
    })

df = pd.DataFrame(records)
output = sys.argv[2]
df.to_csv(output, sep='\t', index=False)
print(f"Created metadata with {len(df)} samples at {output}", file=sys.stderr)
PYTHON "$BATCH_DIR" "$TEMP_METADATA"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate opfi

# Run Sniffles2 stage (6 organisms in parallel)
echo "=========================================="
echo "Running Sniffles2 for IS detection..."
echo "=========================================="

python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
    --metadata "$TEMP_METADATA" \
    --outdir "$BATCH_DIR" \
    --steps sniffles \
    --sniffles-parallel 6 \
    --threads 8

# Run circle detection stage
echo "=========================================="
echo "Running circle detection..."
echo "=========================================="

python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
    --metadata "$TEMP_METADATA" \
    --outdir "$BATCH_DIR" \
    --steps circle \
    --threads 8

# Cleanup
rm -f "$TEMP_METADATA"

echo "=========================================="
echo "Completed: $BATCH_NAME"
echo "Time: $(date)"
echo "=========================================="
