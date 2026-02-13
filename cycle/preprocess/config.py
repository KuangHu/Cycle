"""Default configuration for the preprocessing pipeline."""

# ── Reference genome resolution ──────────────────────────────────────
DEFAULT_REFERENCE_DIR = "data/reference_genomes"

# ── Alignment ────────────────────────────────────────────────────────
DEFAULT_ALIGNMENT_DIR = "data/alignments"
DEFAULT_MINIMAP2_PRESET = "map-ont"
DEFAULT_THREADS = 8
DEFAULT_SORT_MEMORY = "4G"

# ── IS reference ─────────────────────────────────────────────────────
DEFAULT_IS_REFERENCE_DIR = "data/is_reference"
ISFINDER_FASTA_URL = (
    "https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences"
    "/master/IS.fna"
)

# ── FASTQ input ──────────────────────────────────────────────────────
DEFAULT_FASTQ_DIR = "data/sra_downloads"

# ── tldr output ──────────────────────────────────────────────────
DEFAULT_TLDR_OUTPUT_DIR = "data/tldr_output"
