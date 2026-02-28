"""Default configuration for partial circle detection."""

DEFAULT_PARTIAL_CIRCLE_OUTPUT_DIR = "data/partial_circle_output"
DEFAULT_MIN_OVERLAP_EACH_SIDE = 50       # bp aligned on each side of junction
DEFAULT_MIN_CIRCLE_SIZE = 100            # minimum [S, E] span in bp
DEFAULT_MAX_CIRCLE_FRACTION = 0.90       # exclude near-full circles (TH module)
DEFAULT_BREAKPOINT_TOLERANCE = 20        # bp window for merging breakpoints
DEFAULT_MIN_SUPPORTING_READS = 2         # min reads per partial circle call
DEFAULT_MIN_CONSENSUS_LENGTH = 200       # skip short IS elements
DEFAULT_MIN_CONSENSUS_ENTROPY = 1.7      # bits/base; filters low-complexity
