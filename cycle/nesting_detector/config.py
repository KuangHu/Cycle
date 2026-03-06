"""Default configuration for IS-within-IS nesting detection."""

DEFAULT_NESTING_OUTPUT_DIR = "data/nesting_output"
DEFAULT_MM2_PRESET = "asm10"             # minimap2 preset (~1% divergence)
DEFAULT_MIN_IDENTITY = 0.90              # minimum alignment block identity
DEFAULT_MIN_BLOCK_LENGTH = 100           # minimum aligned block (bp)
DEFAULT_MIN_INSERTION_SIZE = 50          # minimum insertion to report (bp)
DEFAULT_FLANKING_PAD = 80               # flanking bp included in extended seq
DEFAULT_MIN_LENGTH_RATIO = 1.02          # host must be >=2% longer than core
