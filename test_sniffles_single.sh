#!/bin/bash
# Test Sniffles2 runner on a single organism using existing batch_000 data
set -e

echo "Testing Sniffles2 on Achromobacter xylosoxidans (1 sample)"

# Create test metadata with just one organism
TEST_META=/tmp/test_achromobacter.tsv
head -1 /home/kuangh/tools/Cycle/data/batches/batch_000.tsv > $TEST_META
grep -i "achromobacter xylosoxidans" /home/kuangh/tools/Cycle/data/batches/batch_000.tsv >> $TEST_META

echo "Created test metadata: $TEST_META"
wc -l $TEST_META

# Use existing batch_000 directories for alignments and reference genomes
# Only create sniffles_output in a new location
BASE_DIR=/groups/rubin/projects/kuang/out/IS_cycle/batch_000
TEST_SNIFFLES_DIR=$BASE_DIR/sniffles_output_test

echo
echo "Running Sniffles2 (using batch_000 BAMs and references)..."
python /home/kuangh/tools/Cycle/scripts/run_pipeline.py \
    --metadata $TEST_META \
    --align-dir $BASE_DIR/alignments \
    --ref-dir $BASE_DIR/reference_genomes \
    --sniffles-dir $TEST_SNIFFLES_DIR \
    --threads 4 \
    --steps sniffles

echo
echo "=== Sniffles2 output ==="
ls -lh $TEST_SNIFFLES_DIR/achromobacter_xylosoxidans/

echo
echo "=== Table preview (first 20 lines) ==="
head -20 $TEST_SNIFFLES_DIR/achromobacter_xylosoxidans/achromobacter_xylosoxidans.table.txt

echo
echo "=== Compare with tldr output ==="
echo "tldr insertions:"
tail -n +2 $BASE_DIR/tldr_output/achromobacter_xylosoxidans/achromobacter_xylosoxidans.table.txt | wc -l
echo "Sniffles2 insertions:"
tail -n +2 $TEST_SNIFFLES_DIR/achromobacter_xylosoxidans/achromobacter_xylosoxidans.table.txt | wc -l

echo
echo "Test complete!"
