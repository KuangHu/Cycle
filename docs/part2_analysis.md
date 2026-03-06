# Part 2: Analysis Modules

Custom analysis scripts that run on outputs from the batch pipeline (Part 1).
These target specific IS families, subsets, or cross-batch questions. Run
interactively or via ad-hoc SLURM jobs — not part of the automated batch
submission.

## Overview

```
Batch pipeline output (Part 1)
  │
  ├─ IS110 Filter ──→ IS110 element set (DEDD+Tnp20 domain requirement)
  │     │
  │     ├─ Nesting Detector ──→ IS-within-IS insertion events
  │     │     └─ Nesting Visualizer
  │     │
  │     └─ Visualization ──→ PNG + GBK per element (domains + partial circles)
  │
  └─ (future modules: IS200, ISAs1, cross-batch statistics, ...)
```

## Modules

### 1. IS110 Filter

Identifies IS elements containing IS110-family transposase. Requires **both**
DEDD nuclease and Tnp20 transposase domains:
- **Case 1:** Single ORF with both DEDD and Tnp20 HMM hits
- **Case 2:** Two adjacent ORFs (within 300bp), one with DEDD and one with Tnp20

Also annotates each ORF in the output records with its domain hits (`"domains": ["DEDD", "Tnp20"]`).

**Script:** `scripts/run_is110_filter.py`
**Module:** `cycle/is110_filter/`

```bash
# Full run (with hmmsearch):
python scripts/run_is110_filter.py \
    --protein-fasta /path/to/system_clustering_batch_NNN/all_proteins.faa \
    --input-dir /path/to/batch_NNN/is_formatter_output \
    --output-dir /path/to/is110_all_batch_NNN \
    --min-tail-head 0 \
    --cpus 8

# Reuse existing hmmsearch results (faster):
python scripts/run_is110_filter.py \
    --protein-fasta UNUSED \
    --input-dir /path/to/batch_NNN/is_formatter_output \
    --output-dir /path/to/is110_all_batch_NNN \
    --min-tail-head 0 \
    --skip-hmmsearch
```

**Key arguments:**
| Flag | Default | Description |
|------|---------|-------------|
| `--protein-fasta` | required | `all_proteins.faa` from system clustering |
| `--input-dir` | required | Formatter output dir with `*_is_records_guide.json` |
| `--output-dir` | required | Where to write filtered records + HMM tblout files |
| `--min-tail-head` | 1 | Set to 0 to include elements without full circle evidence |
| `--evalue` | 1e-5 | hmmsearch E-value threshold |
| `--max-orf-gap` | 300 | Max gap (bp) between ORFs for case 2 adjacency |
| `--skip-hmmsearch` | false | Reuse existing `DEDD_hits.tbl` / `Tnp20_hits.tbl` |

**Output:**
- `is110_circular_records.json` — Filtered records with domain-annotated ORFs
- `is110_circular_summary.tsv` — One-line-per-record summary
- `DEDD_hits.tbl`, `Tnp20_hits.tbl` — hmmsearch tblout files

**HMM paths:** Configured in `cycle/is110_filter/filter.py`:
- DEDD: `/home/kuangh/scripts/IS110/hmm/DEDD.hmm`
- Tnp20: `/home/kuangh/scripts/IS110/hmm/Tnp20.hmm`

---

### 2. Visualization

Generates PNG diagrams and GenBank files for IS elements. Shows:
- Flanking regions (light blue upstream, light yellow downstream)
- ORFs as directional arrows (gray default, **red** DEDD, **orange** Tnp20, **purple** both)
- Noncoding regions (light gray)
- Guide RNA alignment hits (red = upstream flanking, blue = downstream flanking)
- Partial circle regions (**teal** bands with dotted boundary lines)
- Circle evidence in title (TH, GH, TG read counts)

**Scripts:**
- `scripts/run_visualizer.py` — Batch visualization (all samples)
- Direct Python usage for custom subsets (see example below)

**Module:** `cycle/visualizer/`

```python
from cycle.visualizer import ISElementVisualizer, ISElementGenBank

viz = ISElementVisualizer()
gbk = ISElementGenBank()

# With partial circle annotations:
viz.save_element_png(
    element=record,
    alignments=record.get("guide_hits", []),
    output_path="output.png",
    circle_info=record.get("circle_evidence"),
    partial_circles=partial_circle_calls,  # list of dicts or None
)

gbk.save_genbank(
    element=record,
    alignments=record.get("guide_hits", []),
    output_path="output.gbk",
    circle_info=record.get("circle_evidence"),
    partial_circles=partial_circle_calls,
)
```

**Partial circle dict format** (from `_partial_circle_summary.json`):
```json
{
    "circle_start": 690,     // 0-based position on IS element
    "circle_end": 2045,
    "circle_size": 1355,
    "circle_fraction": 0.31,
    "n_supporting_reads": 37
}
```

**Domain colors in PNG:**
| Color | Meaning |
|-------|---------|
| Gray `#bbbbbb` | ORF without domain annotation |
| Red `#e41a1c` | ORF with DEDD domain only |
| Orange `#ff7f00` | ORF with Tnp20 domain only |
| Purple `#984ea3` | ORF with both DEDD + Tnp20 |
| Teal `#26a69a` | Partial circle region |

---

### 3. Nesting Detector

Detects IS-within-IS insertion events by pairwise minimap2 alignment. A "host"
element contains an insertion relative to a shorter "core" element.

**Script:** `scripts/run_nesting_detector.py`
**Module:** `cycle/nesting_detector/`

```bash
python scripts/run_nesting_detector.py \
    --records /path/to/is110_circular_records.json \
    --output-dir /path/to/nesting_output \
    --min-identity 0.90 \
    --min-insertion-size 50 \
    --flanking-pad 80 \
    --threads 8
```

**Key arguments:**
| Flag | Default | Description |
|------|---------|-------------|
| `--records` | required | IS records JSON (e.g., IS110 filtered output) |
| `--output-dir` | required | Output directory |
| `--min-identity` | 0.90 | Minimum alignment identity |
| `--min-insertion-size` | 50 | Minimum insertion size (bp) |
| `--min-block-length` | 100 | Minimum aligned block length |
| `--flanking-pad` | 80 | Flanking bp added to each side for alignment |
| `--min-length-ratio` | 1.02 | Host must be ≥ this × core length |
| `--threads` | 8 | minimap2 threads |

**Algorithm:**
1. Build extended sequences (upstream_flank + IS + downstream_flank)
2. All-vs-all minimap2 alignment (`-x asm10 -X -c --eqx`)
3. Decompose CIGAR at large indels into alignment blocks
4. Detect collinear block arrangements with insertions in the host

**Output:**
- `nesting_events.json` — Detailed event records with aligned blocks
- `nesting_events.tsv` — Flat summary

---

### 4. Nesting Visualizer

Draws paired diagrams for nesting events (host + core) with aligned blocks
(green) and insertion regions (coral).

**Script:** `scripts/visualize_nesting_events.py`

```bash
python scripts/visualize_nesting_events.py \
    --events /path/to/nesting_events.json \
    --records /path/to/is110_circular_records.json \
    --output-dir /path/to/nesting_vis \
    --max-events 50
```

---

## Typical Workflow

### IS110 analysis on a new batch

```bash
BATCH=batch_021
BASE=/groups/rubin/projects/kuang/out/IS_cycle

# 1. Run IS110 filter (requires step 5 clustering output for proteins)
python scripts/run_is110_filter.py \
    --protein-fasta $BASE/system_clustering_${BATCH}/all_proteins.faa \
    --input-dir $BASE/$BATCH/is_formatter_output \
    --output-dir $BASE/is110_all_${BATCH} \
    --min-tail-head 0

# 2. Cross-reference with partial circle data
python3 -c "
import json, os
from collections import defaultdict

with open('$BASE/is110_all_${BATCH}/is110_circular_records.json') as f:
    records = json.load(f)

pc_dir = '$BASE/$BATCH/partial_circle_output'
pc_by_is = defaultdict(list)
for r in records:
    sid = r['sample_id']
    sj = os.path.join(pc_dir, sid, f'{sid}_partial_circle_summary.json')
    if not os.path.exists(sj): continue
    with open(sj) as f:
        for c in json.load(f):
            if c['is_id'] == r['is_id']:
                pc_by_is[c['is_id']].append(c)

print(f'IS110 elements: {len(records)}')
print(f'With partial circles: {len(pc_by_is)}')
"

# 3. Visualize all IS110 elements
python3 -c "
import sys, json, os
from collections import defaultdict
sys.path.insert(0, '.')
from cycle.visualizer import ISElementVisualizer, ISElementGenBank

with open('$BASE/is110_all_${BATCH}/is110_circular_records.json') as f:
    records = json.load(f)

pc_dir = '$BASE/$BATCH/partial_circle_output'
pc_by_is = defaultdict(list)
for r in records:
    sid = r['sample_id']
    sj = os.path.join(pc_dir, sid, f'{sid}_partial_circle_summary.json')
    if not os.path.exists(sj): continue
    seen = set()
    with open(sj) as f:
        for c in json.load(f):
            if c['is_id'] in {r['is_id'] for r in records}:
                key = (c['circle_start'], c['circle_end'])
                if key not in seen:
                    seen.add(key)
                    pc_by_is[c['is_id']].append(c)

outdir = '$BASE/is110_all_${BATCH}/vis_all'
os.makedirs(outdir, exist_ok=True)
viz = ISElementVisualizer()
gbk = ISElementGenBank()

for rec in records:
    is_id = rec['is_id']
    sid = rec['sample_id']
    pcs = pc_by_is.get(is_id) or None
    tag = f'{sid}_{is_id[:12]}'
    viz.save_element_png(rec, rec.get('guide_hits', []),
        os.path.join(outdir, f'{tag}.png'),
        circle_info=rec.get('circle_evidence'), partial_circles=pcs)
    gbk.save_genbank(rec, rec.get('guide_hits', []),
        os.path.join(outdir, f'{tag}.gbk'),
        circle_info=rec.get('circle_evidence'), partial_circles=pcs)
print(f'Drew {len(records)} elements')
"

# 4. Run nesting detector
python scripts/run_nesting_detector.py \
    --records $BASE/is110_all_${BATCH}/is110_circular_records.json \
    --output-dir $BASE/is110_all_${BATCH}/nesting
```

## Data Locations

All outputs live under `/groups/rubin/projects/kuang/out/IS_cycle/`:

```
IS_cycle/
  batch_000/                        # Per-batch pipeline output (Part 1)
    sra_downloads/
    reference_genomes/
    alignments/
    sniffles_output/
    circle_output/
    is_formatter_output/            # Final per-sample records
    partial_circle_output/          # Per-sample partial circle calls
  ...
  batch_094/

  system_clustering_batch_000/      # Per-batch clustering (Part 1 step 5)
  novelty_batch_000/                # Per-batch novelty scores (Part 1 step 5)
  ...

  is110_all_batch_000/              # IS110 analysis (Part 2)
    is110_circular_records.json
    DEDD_hits.tbl
    Tnp20_hits.tbl
    vis_all/                        # PNG + GBK visualizations
    nesting/                        # Nesting detection output
  ...
```
