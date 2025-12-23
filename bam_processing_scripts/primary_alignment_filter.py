#!/usr/bin/env python3
"""
primary_alignment_filter.py

Split a SAM file into:
  - Primary C. elegans alignments
  - Non-primary alignments (secondary / supplementary)

The FLAG field is used to determine if an alignment is primary:
- Primary: neither 0x100 (secondary) nor 0x800 (supplementary) bits are set.

Usage
-----
    python primary_alignment_filter.py <input.sam>

Requirements
------------
- SAM file with a valid integer FLAG column.
- Input filename is assumed to end with "_celegans.sam"; output names are:
      *_cel.sam         (primary alignments)
      *_not_primary.sam (non-primary)
"""

import sys


def process_sam_file(sam_file_path):
    """
    Process a SAM file and split it into:
    - C. elegans primary alignments
    - Non-primary alignments (secondary or supplementary)

    The FLAG field is used to determine if an alignment is primary.
    """

    # Prepare output filenames based on the input SAM file name
    base_name = sam_file_path.replace("_celegans.sam", "")
    celegans_file = f"{base_name}_cel.sam"
    not_primary_file = f"{base_name}_not_primary.sam"

    # Open output files for writing results
    with open(sam_file_path, 'r') as sam_file, \
            open(celegans_file, 'w') as celegans_out, \
            open(not_primary_file, 'w') as not_primary_out:

        for line in sam_file:
            # Write headers directly to all output files for BAM validity
            if line.startswith('@'):
                celegans_out.write(line)
                not_primary_out.write(line)
                continue

            # Split the line into fields to extract FLAG
            fields = line.strip().split("\t")
            try:
                flag = int(fields[1])  # Convert FLAG to integer for bitmask comparison
            except ValueError:
                print(f"Invalid FLAG value in line: {line}")
                continue

            # Check if the alignment is primary using bitwise operations
            is_primary = not (flag & 256 or flag & 2048)

            # Write the line to the appropriate output file
            if is_primary:
                celegans_out.write(line)
            else:
                not_primary_out.write(line)

    print(f"Processing complete for {sam_file_path}")
    print(f"Primary alignments saved to: {celegans_file}")
    print(f"Non-primary alignments saved to: {not_primary_file}")


# Command-line execution
if __name__ == "__main__":
    # Ensure the script is provided with a valid SAM file
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input_sam_file>")
        sys.exit(1)

    # Process the provided SAM file
    input_sam_file = sys.argv[1]
    process_sam_file(input_sam_file)
