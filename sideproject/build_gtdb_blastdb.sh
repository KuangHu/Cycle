#!/bin/bash
#SBATCH --job-name=gtdb_blastdb
#SBATCH --output=/home/kuangh/logs/gtdb_blastdb_%j.out
#SBATCH --error=/home/kuangh/logs/gtdb_blastdb_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --qos=standard
#SBATCH --partition=standard
#SBATCH --exclude=node-48-256g-2,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-256g-12,node-48-256g-13,node-48-256g-14,node-48-256g-17,node-48-256g-18,node-48-256g-19,node-48-256g-20

DB_DIR=/groups/rubin/projects/kuang/db/gtdb_blastdb

echo "=== Building BLAST database (combined FASTA already exists) ==="
ls -lh $DB_DIR/gtdb_reps_r220.fna

# Clean up chunks from previous run
rm -f $DB_DIR/chunk_* $DB_DIR/genome_list.txt

echo "=== Running makeblastdb ==="
makeblastdb \
    -in $DB_DIR/gtdb_reps_r220.fna \
    -dbtype nucl \
    -out $DB_DIR/gtdb_reps_r220 \
    -title "GTDB r220 representative genomes" \
    -parse_seqids

echo "=== Done ==="
ls -lh $DB_DIR/gtdb_reps_r220.*
