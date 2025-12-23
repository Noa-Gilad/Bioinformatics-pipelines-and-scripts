#!/usr/bin/env python3
"""
ambi_chr_sep.py

Split C. elegans SAM alignments into X-chromosome vs autosome/ambiguous.

This script expects a SAM file in which all records for the same read
(QNAME) are contiguous (i.e., sorted by QNAME). For each read it checks
the reference name (RNAME) across all its alignments and classifies it as:

  * X-only        – all alignments on chromosome X
  * Ambiguous_chr – at least one alignment on X and at least one on
                    a non-X chromosome (autosome)

Two SAM files are written, derived from the input basename:
  <basename>_celegans_chr.sam
  <basename>_cel_ambiguous_chr.sam

Usage
-----
    python ambi_chr_sep.py <input_celegans.sam>

Notes
-----
- SAM headers (lines starting with '@') are copied to both outputs.
- The exact chromosome naming convention ("X" vs "chrX", etc.) should
  be adjusted in the code if your reference uses different names.
"""

import sys


def process_sam_file(sam_file_path):
    """
    Process a SAM file (already sorted by QNAME) and split it into:
    - C. elegans X chromosome only alignments
    - Ambiguous reads (X + autosome)
    """
    # Prepare output filenames based on the input SAM file name
    base_name = sam_file_path.replace("_celegans.sam", "")
    celegans_file = f"{base_name}_celegans_chr.sam"
    ambiguous_file = f"{base_name}_cel_ambiguous_chr.sam"

    # Open output files for writing results
    with open(sam_file_path, 'r') as sam_file, \
            open(celegans_file, 'w') as celegans_out, \
            open(ambiguous_file, 'w') as ambiguous_out:

        # Initialize variables for processing
        chunk = []           # Store current QNAME block
        current_qname = None  # Track the current read group name

        for line in sam_file:
            # Write headers directly to all output files for BAM validity
            if line.startswith('@'):
                celegans_out.write(line)
                ambiguous_out.write(line)
                continue

            # Split the line into fields to extract QNAME and RNAME
            fields = line.strip().split("\t")
            qname = fields[0]  # QNAME (Read Name)

            # If the QNAME changes, process the previous chunk
            if qname != current_qname and chunk:
                # Classify the previous read group based on RNAME
                X_present = any("X" in l.split("\t")[2] for l in chunk)
                autosomes_present = any("X" not in l.split("\t")[2] for l in chunk)

                # Write to the appropriate file based on classification
                if X_present and autosomes_present:
                    ambiguous_out.writelines(chunk)
                else:
                    celegans_out.writelines(chunk)

                # Reset the chunk for the next QNAME
                chunk = []

            # Add the current line to the chunk
            chunk.append(line)
            current_qname = qname

        # Process the final chunk (last read group)
        if chunk:
            X_present = any("X" in l.split("\t")[2] for l in chunk)
            autosomes_present = any("X" not in l.split("\t")[2] for l in chunk)

            if X_present and autosomes_present:
                ambiguous_out.writelines(chunk)
            else:
                celegans_out.writelines(chunk)

    print(f"Processing complete for {sam_file_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python ambi_chr_sep.py <input_sam_file>")
        sys.exit(1)

    input_sam_file = sys.argv[1]
    process_sam_file(input_sam_file)
