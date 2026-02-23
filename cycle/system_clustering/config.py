"""Default parameters for system clustering."""

# MMseqs2 clustering
DEFAULT_MIN_SEQ_ID = 0.3
DEFAULT_COVERAGE = 0.8
DEFAULT_COV_MODE = 0  # bidirectional coverage
DEFAULT_MMSEQS_THREADS = 8

# Community detection (Louvain)
DEFAULT_LOUVAIN_RESOLUTION = 1.0

# Jaccard similarity threshold for graph edges
DEFAULT_JACCARD_SIM_THRESHOLD = 0.3

# Variant analysis
DEFAULT_FLANKING_EDIT_THRESHOLD = 10  # max edit distance for same L1 group
