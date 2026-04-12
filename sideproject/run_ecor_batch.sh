#!/bin/bash
#SBATCH --job-name=ecor_batch
#SBATCH --output=/home/kuangh/logs/ecor_batch_%j.out
#SBATCH --error=/home/kuangh/logs/ecor_batch_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --qos=standard
#SBATCH --partition=standard

cd /home/kuangh/tools/Cycle
conda run -n opfi python sideproject/run_ecor_batch.py \
    --fastq-dir /groups/rubin/projects/kuang/out/IS110/ECOR_batch/fastq \
    --ref-dir /groups/rubin/projects/kuang/out/IS110/ECOR_batch/reference_genomes \
    --output-dir /groups/rubin/projects/kuang/out/IS110/ECOR_batch/output \
    --threads 48
