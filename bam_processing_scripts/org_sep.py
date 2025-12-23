#!/usr/bin/env python3
"""
org_sep.py

Split a SAM file (sorted by QNAME) into three outputs:
  - One file containing reads aligned to the E. coli genome
  - One file containing reads aligned to the C. elegans genome
  - One file containing reads that map to both (ambiguous)

Usage
-----
    python org_sep.py <input.sam>

Requirements
------------
- Input SAM must be sorted by QNAME (all records of the same read are adjacent).
- SAM records for E. coli use an RNAME that contains the string "Chromosome".
  C. elegans records are assumed to have RNAME values that do *not* contain
  "Chromosome".

Notes
-----
- Output filenames are derived from the input name by stripping a specific
  pipeline suffix and appending:
      *_ecoli.sam, *_celegans.sam, *_ambiguous.sam
"""

import sys


def process_sam_file(sam_file_path):
    """
    Process a SAM file (already sorted by QNAME) and split it into:
    - E. coli (RNAME == "Chromosome")
    - C. elegans (RNAME != "Chromosome")
    - Ambiguous (both references present in a read group)
    """

    # Prepare output filenames based on the input SAM file name
    base_name = sam_file_path.replace(
        "_Caenorhabditis_elegans_miRNA-Seq_outputAligned.sortedByCoord.out.sam",
        ""
    )
    ecoli_file = f"{base_name}_ecoli.sam"
    celegans_file = f"{base_name}_celegans.sam"
    ambiguous_file = f"{base_name}_ambiguous.sam"

    # Open output files for writing results
    with open(sam_file_path, 'r') as sam_file, \
            open(ecoli_file, 'w') as ecoli_out, \
            open(celegans_file, 'w') as celegans_out, \
            open(ambiguous_file, 'w') as ambiguous_out:

        # Initialize variables for processing
        chunk = []  # Store current QNAME block
        current_qname = None  # Track the current read group name

        for line in sam_file:
            # Write headers directly to all output files for BAM validity
            if line.startswith('@'):
                ecoli_out.write(line)
                celegans_out.write(line)
                ambiguous_out.write(line)
                continue

            # Split the line into fields to extract QNAME and RNAME
            fields = line.strip().split("\t")
            qname = fields[0]  # QNAME (Read Name)

            # If the QNAME changes, process the previous chunk
            if qname != current_qname and chunk:
                # Classify the previous read group based on RNAME
                ecoli_present = any(
                    "Chromosome" in ln.split("\t")[2] for ln in chunk
                )
                celegans_present = any(
                    "Chromosome" not in ln.split("\t")[2] for ln in chunk
                )

                # Write to the appropriate file based on classification
                if ecoli_present and not celegans_present:
                    ecoli_out.writelines(chunk)
                elif celegans_present and not ecoli_present:
                    celegans_out.writelines(chunk)
                else:
                    ambiguous_out.writelines(chunk)

                # Reset the chunk for the next QNAME
                chunk = []

            # Add the current line to the chunk
            chunk.append(line)
            current_qname = qname

        # Process the final chunk (last read group)
        if chunk:
            ecoli_present = any(
                "Chromosome" in ln.split("\t")[2] for ln in chunk
            )
            celegans_present = any(
                "Chromosome" not in ln.split("\t")[2] for ln in chunk
            )

            if ecoli_present and not celegans_present:
                ecoli_out.writelines(chunk)
            elif celegans_present and not ecoli_present:
                celegans_out.writelines(chunk)
            else:
                ambiguous_out.writelines(chunk)

    print(f"Processing complete for {sam_file_path}")


if __name__ == "__main__":
    # Ensure the script is provided with a valid SAM file
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_sam_file>")
        sys.exit(1)

    # Process the provided SAM file
    input_sam_file = sys.argv[1]
    process_sam_file(input_sam_file)
