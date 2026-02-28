# IS_cycle Data Structure Guide

Root: `/groups/rubin/projects/kuang/out/IS_cycle/`

## Per-Batch Directories

Each `batch_NNN/` contains per-sample pipeline outputs only. All data within a
batch directory is produced by processing individual samples independently.

```
batch_NNN/
├── sra_downloads/                 # Raw FASTQs from SRA
│   ├── SRR12345678.fastq.gz
│   ├── ERR12345678_1.fastq.gz
│   └── ...
│
├── reference_genomes/             # NCBI reference genomes
│   ├── resolve_status.tsv         # Organism → accession mapping
│   ├── GCF_000005845.2/           # One dir per assembly
│   │   ├── *.fna.gz               # Genome FASTA
│   │   └── *.mmi                  # minimap2 index
│   └── ...
│
├── alignments/                    # minimap2 + samtools sort
│   ├── alignment_status.tsv       # Per-sample alignment status
│   ├── SRR12345678.sorted.bam
│   ├── SRR12345678.sorted.bam.bai
│   └── ...
│
├── is_reference/                  # ISfinder reference
│   └── ISfinder_raw.fna           # ISfinder sequences (may be absent
│                                  #   in sniffles-only batches; use
│                                  #   batch_000's copy for novelty)
│
├── sniffles_output/               # Sniffles2 structural variant calls
│   ├── organism_name/             # Keyed by organism (slug)
│   │   ├── organism.table.txt     # IS insertion table
│   │   └── SRR12345678.vcf        # Sniffles VCF
│   └── ...
│
├── circle_output/                 # Full-circle (tail-head) detection
│   ├── organism_name/             # Keyed by organism (slug)
│   │   ├── *_bait.fa              # [IS][IS] concatemer bait
│   │   ├── *_circle_reads.tsv     # Junction-spanning reads
│   │   └── *_circle_summary.tsv   # Per-IS TH read counts
│   └── ...
│
├── is_formatter_output/           # ★ MAIN RESULTS — unified per-sample IS records
│   ├── SRR12345678/               # Keyed by sample accession
│   │   ├── assembly/              # Reconstructed IS element FASTAs
│   │   ├── SRR12345678_is_records.json       # Base IS records
│   │   ├── SRR12345678_is_records.tsv        # Flat TSV version
│   │   ├── SRR12345678_is_records_annotated.json  # + ORF annotation
│   │   └── SRR12345678_is_records_guide.json      # + guide alignment hits
│   └── ...
│
├── partial_circle_output/         # Partial circle (sub-region) detection
│   ├── SRR12345678/               # Keyed by sample accession
│   │   ├── *_is_ref.fa            # Single-copy IS reference
│   │   ├── *_partial_circle_reads.tsv     # Supporting split reads
│   │   ├── *_partial_circle_summary.json  # Clustered calls (detailed)
│   │   ├── *_partial_circle_summary.tsv   # Clustered calls (flat)
│   │   └── *.partial.nsorted.bam          # Name-sorted BAM (intermediate)
│   └── ...
│
├── metadata_for_sniffles.tsv      # Sample metadata for pipeline
└── partial_circle_manifest.tsv    # Sample manifest for partial circle jobs
```

### Key per-sample files

The most important output for downstream analysis is `is_formatter_output/`.
Each sample's `*_is_records_guide.json` is the richest record, containing:

- IS element sequence + coordinates + family classification
- Flanking region sequences (upstream/downstream, typically 80bp)
- ORF annotations (start, end, strand, protein sequence)
- Noncoding region annotations
- Guide alignment hits (flanking-to-noncoding matches)
- Circle evidence (tail-head junction read counts)

## Cross-Batch Analysis Directories

These live at the `IS_cycle/` root level. They can operate on one batch or
combine multiple batches.

```
IS_cycle/
├── system_clustering_batch_NNN/       # Protein clustering (single batch)
│   ├── all_proteins.faa               # All ORF protein sequences
│   ├── mmseqs_clusters/               # MMseqs2 raw clustering output
│   ├── system_clusters.json           # Louvain community clusters
│   ├── system_clusters_summary.tsv    # Per-IS cluster assignments
│   ├── family_ids.json                # Cluster → family mapping
│   ├── transposon_ids.json            # IS → cluster mapping
│   └── incidence_matrix.npz           # Sparse co-occurrence matrix
│
├── system_clustering_batch_000_to_004/  # Combined clustering (example)
│   └── ...                              # Same structure as above
│
├── novelty_batch_NNN/                 # ISfinder novelty annotation
│   ├── query_is_elements.fna          # IS nucleotide sequences for BLAST
│   ├── isfinder_blastdb.*             # Local ISfinder BLAST database
│   ├── blast_results.tsv              # Raw BLAST hits
│   ├── annotated_is_records.json      # Per-IS novelty annotations
│   ├── cluster_novelty.json           # Per-cluster novelty (detailed)
│   ├── cluster_novelty_summary.tsv    # Per-IS novelty scores (flat)
│   └── novel_candidates/              # Curated candidate PNGs + GBKs
│
├── is110_circular_batch_NNN/          # IS110-specific circular analysis
│   ├── is110_circular_records.json
│   ├── is110_circular_summary.tsv
│   ├── DEDD_hits.tbl                  # DEDD nuclease domain hits
│   ├── Tnp20_hits.tbl                 # Tnp20 transposase domain hits
│   └── visualizations/
│
└── ...
```

### Naming convention

- Per-batch: `{analysis}_batch_NNN/` (e.g., `novelty_batch_000/`)
- Combined: `{analysis}_batch_NNN_to_MMM/` (e.g., `system_clustering_batch_000_to_004/`)

## Pipeline Flow

```
sra_downloads → reference_genomes → alignments
                                        │
                              ┌─────────┴─────────┐
                              ▼                     ▼
                      sniffles_output         (tldr_output)
                              │                     │
                              ▼                     ▼
                        circle_output         (legacy path)
                              │
                              ▼
                    is_formatter_output    ← main per-sample results
                       │            │
                       ▼            ▼
              partial_circle    system_clustering  ← cross-batch
                                       │
                                       ▼
                                    novelty        ← cross-batch
```

## Notes

- **Organism-keyed vs sample-keyed**: `sniffles_output/` and `circle_output/`
  are keyed by organism slug (one organism can have multiple samples).
  `is_formatter_output/` and `partial_circle_output/` are keyed by sample
  accession (SRR/ERR/DRR).

- **ISfinder FASTA**: Batches that used the sniffles pipeline (005+) do not
  have `is_reference/`. Use `batch_000/is_reference/ISfinder_raw.fna` when
  running novelty annotation on those batches (pass `--isfinder-fasta`).

- **Alignments are large**: BAM files dominate storage (~700GB per batch).
  They are only needed if re-running sniffles/circle detection. Once
  `is_formatter_output/` is finalized, BAMs can be deleted to save space.

- **Partial circle BAMs**: Each sample's `*.partial.nsorted.bam` in
  `partial_circle_output/` can also be deleted after the summary files are
  produced, to save space.
