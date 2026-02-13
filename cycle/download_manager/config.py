"""Default search queries and configuration for SRA downloads."""

# NCBI Entrez base URL
ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Default SRA search query: bacterial nanopore WGS genomic reads
DEFAULT_SRA_QUERY = (
    '(("Bacteria"[Organism] OR bacteria[All Fields])'
    ' AND "Oxford Nanopore"[Platform]'
    ' AND "WGS"[Strategy]'
    ' AND "genomic"[Source])'
)

# Entrez DB
SRA_DB = "sra"

# Max records per Entrez request (NCBI hard limit is 10000)
BATCH_SIZE = 500

# Seconds between Entrez requests (NCBI asks for <=3 req/s without API key)
REQUEST_DELAY = 0.4

# Kingfisher download settings
DEFAULT_DOWNLOAD_DIR = "data/sra_downloads"
DEFAULT_DOWNLOAD_METHODS = ["ena-ftp", "aws-http", "prefetch"]
DEFAULT_OUTPUT_FORMAT = "fastq.gz"
DEFAULT_DOWNLOAD_THREADS = 8
DEFAULT_EXTRACTION_THREADS = 8
