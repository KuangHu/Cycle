#!/bin/bash
#SBATCH --job-name=is110_blastn
#SBATCH --output=/home/kuangh/logs/is110_blastn_%j.out
#SBATCH --error=/home/kuangh/logs/is110_blastn_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --qos=standard
#SBATCH --partition=standard

DB=/groups/rubin/projects/kuang/db/gtdb_blastdb/gtdb_reps_r220
QUERY=/groups/rubin/projects/kuang/out/IS_cycle/is110_all_006_026/is110_consensus.fna
OUT_DIR=/groups/rubin/projects/kuang/out/IS_cycle/is110_split_gtdb
mkdir -p $OUT_DIR

echo "=== BLASTn: 2276 IS110 vs GTDB ==="
blastn \
    -query $QUERY \
    -db $DB \
    -out $OUT_DIR/is110_vs_gtdb.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-10 \
    -num_threads 48 \
    -max_target_seqs 500 \
    -max_hsps 5

echo "=== Results ==="
echo "Hits: $(wc -l < $OUT_DIR/is110_vs_gtdb.tsv)"
echo "Unique IS110 with hits: $(cut -f1 $OUT_DIR/is110_vs_gtdb.tsv | sort -u | wc -l)"
echo "Unique genomes hit: $(cut -f2 $OUT_DIR/is110_vs_gtdb.tsv | sort -u | wc -l)"

echo "=== Finding split hits ==="
conda run -n opfi python /home/kuangh/tools/Cycle/sideproject/find_splits_from_blast.py \
    --blast-results $OUT_DIR/is110_vs_gtdb.tsv \
    --output-dir $OUT_DIR

echo "=== Done ==="
