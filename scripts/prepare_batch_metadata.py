#!/usr/bin/env python3
"""Generate metadata for Sniffles2 + circle detection from alignment status."""

import pandas as pd
import sys
from pathlib import Path

def get_organism_from_reference(ref_dir):
    """Extract organism name from reference genome directory."""
    # Look for the genomic.fna file and try to parse organism from it
    fasta_files = list(ref_dir.glob("*.fna"))
    if not fasta_files:
        return ref_dir.name  # Use accession as fallback

    # Read first line of FASTA to get organism name
    try:
        with open(fasta_files[0]) as f:
            header = f.readline().strip()
            # Parse organism from header like: >NC_123456.1 Escherichia coli strain...
            if ' ' in header:
                parts = header.split('[')
                if len(parts) > 1:
                    org_part = parts[1].split(']')[0]
                    if org_part:
                        return org_part
                # Try simpler parsing
                words = header[1:].split()
                if len(words) >= 2:
                    # Take first two words as organism (genus species)
                    return f"{words[0]} {words[1]}"
    except:
        pass

    return ref_dir.name

def main(batch_dir):
    batch_path = Path(batch_dir)
    alignment_status = batch_path / "alignments" / "alignment_status.tsv"

    if not alignment_status.exists():
        print(f"Error: {alignment_status} not found", file=sys.stderr)
        sys.exit(1)

    # Read alignment status
    df = pd.read_csv(alignment_status, sep='\t')

    # Extract organism from reference path
    organisms = []
    for ref_path in df['reference']:
        ref_path = Path(ref_path)
        ref_dir = ref_path.parent
        organism = get_organism_from_reference(ref_dir)
        organisms.append(organism)

    # Create metadata
    metadata = pd.DataFrame({
        'srr_accession': df['sample_id'],
        'organism': organisms,
        'reference': df['reference'],
    })

    # Write metadata
    output = batch_path / "metadata_for_sniffles.tsv"
    metadata.to_csv(output, sep='\t', index=False)

    print(f"Created metadata at {output}")
    print(f"  Samples: {len(metadata)}")
    print(f"  Organisms: {metadata['organism'].nunique()}")

    # Show sample
    print("\nFirst few rows:")
    print(metadata.head())

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: prepare_batch_metadata.py <batch_dir>")
        sys.exit(1)
    main(sys.argv[1])
