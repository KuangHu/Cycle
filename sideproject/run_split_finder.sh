#!/bin/bash
#SBATCH --job-name=is110_split
#SBATCH --output=/home/kuangh/logs/is110_split_%j.out
#SBATCH --error=/home/kuangh/logs/is110_split_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --qos=standard
#SBATCH --partition=standard

cd /home/kuangh/tools/Cycle
conda run -n opfi python sideproject/is110_split_finder.py \
    --is110-fasta /groups/rubin/projects/kuang/out/IS_cycle/is110_all_006_026/is110_consensus.fna \
    --genome-dir /groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether \
    --output-dir /groups/rubin/projects/kuang/out/IS_cycle/is110_split_gtdb \
    --threads 48 \
    --batch-size 2000
