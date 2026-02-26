"""Default parameters for ISfinder novelty annotation."""

# ── BLAST search ─────────────────────────────────────────────────────
DEFAULT_BLAST_EVALUE = 1e-5
DEFAULT_BLAST_MAX_TARGET_SEQS = 5
DEFAULT_BLAST_THREADS = 8

# ── Novelty scoring weights ─────────────────────────────────────────
WEIGHT_DIVERGENCE = 0.5
WEIGHT_DIVERSITY = 0.3
WEIGHT_MOSAIC = 0.2

# ── Divergence thresholds ───────────────────────────────────────────
# Below this mean %identity the divergence score saturates at 1.0
PIDENT_FLOOR = 50.0
# At or above this %identity the divergence score is 0.0
PIDENT_CEILING = 95.0

# ── Novelty class boundaries ────────────────────────────────────────
NOVEL_THRESHOLD = 0.7
DIVERGENT_THRESHOLD = 0.4

# ── Minimum query coverage to count a BLAST hit ─────────────────────
MIN_QUERY_COVERAGE = 0.3

# ── Output filenames ────────────────────────────────────────────────
ANNOTATED_RECORDS_JSON = "annotated_is_records.json"
CLUSTER_NOVELTY_TSV = "cluster_novelty_summary.tsv"
CLUSTER_NOVELTY_JSON = "cluster_novelty.json"
