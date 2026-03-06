# Setting up Cycle on Lawrence Cluster

This guide covers what needs to change to run the Cycle pipeline on Lawrence (LRC) instead of Savio.

## 1. Clone and install

```bash
git clone <repo-url> ~/tools/Cycle
cd ~/tools/Cycle

# Create conda env with all dependencies
conda create -n opfi python=3.11
conda activate opfi
pip install pandas biopython pysam
conda install -c bioconda minimap2 samtools kingfisher sniffles=2.7.2 prodigal hmmer mafft exonerate mmseqs2 blast
```

## 2. Files that need cluster-specific edits

### submit_one_batch.sh (main submission script)

Lines ~70-77 — change these variables:
```bash
CYCLE_DIR=$HOME/tools/Cycle                    # was /home/kuangh/tools/Cycle
OUTROOT=$HOME/scratch/IS_cycle                 # was /groups/rubin/projects/kuang/out/IS_cycle
CONDA_ENV=opfi
ENV_SETUP="export PATH=$HOME/miniconda3/envs/${CONDA_ENV}/bin:\$PATH"  # adjust conda path
```

Lines with `--partition=standard --qos=standard` — change to Lawrence partition/QOS:
```bash
--partition=lr_normal --qos=lr_normal
# or whatever partition you have access to on Lawrence
```

Line ~75 — ISfinder reference path:
```bash
ISFINDER=$OUTROOT/batch_000/is_reference/ISfinder_raw.fna
# This gets created by the pipeline. For the first batch, it downloads automatically.
```

### submit_all.sh (old submission script, same changes)

Same variables at lines ~60-73 and partition/QOS settings.

### cycle/is110_filter/filter.py (lines 28-29)

Hardcoded HMM paths for IS110 domain detection:
```python
DEFAULT_DEDD_HMM = "/home/kuangh/scripts/IS110/hmm/DEDD.hmm"
DEFAULT_TNP20_HMM = "/home/kuangh/scripts/IS110/hmm/Tnp20.hmm"
```

Options:
- Copy the HMM files into the repo at `data/hmm/` and update these defaults
- Or pass `--dedd-hmm` / `--tnp20-hmm` on the CLI when running IS110 filter

### scripts/check_status.py (lines 147, 151)

Default paths in argparse:
```python
"--outroot", default="/groups/rubin/projects/kuang/out/IS_cycle"
"--batchdir", default="/home/kuangh/tools/Cycle/data/batches"
```

Change to Lawrence paths, or always pass them explicitly.

## 3. Files that are already portable (no changes needed)

- All Python modules in `cycle/` — paths come from CLI args
- `cycle/preprocess/config.py` — uses relative paths
- `scripts/run_pipeline.py` — everything via `--outdir`, `--metadata`
- `scripts/run_orf_annotator.py`, `run_guide_finder.py`, etc. — all CLI-driven
- Batch metadata TSVs in `data/batches/` — just sample accessions, no paths

## 4. Data directory structure

Create the output root and it will be populated by the pipeline:
```bash
mkdir -p $HOME/scratch/IS_cycle
```

The pipeline creates per-batch directories:
```
IS_cycle/
  batch_000/
    sra_downloads/       # FASTQs (deleted after processing)
    reference_genomes/   # one per organism (kept)
    is_formatter_output/ # final IS records + JSONs (kept)
    chunks/              # per-node working dirs (deleted by cleanup phase)
  system_clustering_batch_000/
  novelty_batch_000/
```

## 5. SLURM differences

| Setting | Savio | Lawrence |
|---------|-------|----------|
| Partition | `standard` | `lr_normal` (check with `sinfo`) |
| QOS | `standard` | `lr_normal` (check with `sacctmgr show qos`) |
| Max submitted | 200 | check QOS limits |
| Max running | 10 | check QOS limits |
| CPUs per node | 48 | check with `sinfo -N -l` |
| Memory per node | 192G | check node specs |

Adjust `--cpus-per-task` and `--mem` in submit_one_batch.sh if nodes are different sizes.

## 6. Quick test

Run a small test to verify everything works:
```bash
# Create a test batch with 2-3 samples
head -4 data/batches/batch_006.tsv > data/batches/batch_test.tsv

# Run pipeline manually (not via SLURM) on one sample
conda activate opfi
python scripts/run_pipeline.py \
  --metadata data/batches/batch_test.tsv \
  --outdir $HOME/scratch/IS_cycle/batch_test \
  --steps download resolve index align sniffles circle format \
  --parallel 1 --threads 8
```

## Summary of changes needed

1. **5 path variables** in `submit_one_batch.sh` (CYCLE_DIR, OUTROOT, conda, partition, QOS)
2. **2 HMM paths** in `cycle/is110_filter/filter.py` (or use CLI flags)
3. **2 default paths** in `scripts/check_status.py` (or pass explicitly)
4. **Resource limits** in sbatch calls if node sizes differ
