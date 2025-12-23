#!/usr/bin/env python3
import os

"""
samtools_batch_utils.py

Description:
    Convenience functions for running common samtools operations on BAM files
    in a directory:
        - BAM → SAM
        - SAM → BAM
        - sorting BAM
        - indexing BAM
        - flagstat statistics

By default, the script:
    - Loads samtools via "module load samtools"
    - Collects all *.out.bam files under BASE_STAR_OUT_DIR

Usage:
    1. Edit BASE_STAR_OUT_DIR below to point to your STAR output directory.
    2. Run this script to loop over BAMs and apply the desired functions
       (the bottom loop currently runs BamToSam on each BAM).

Requirements:
    - samtools installed and available in PATH or via "module load samtools".

Paths requiring user configuration:
    - BASE_STAR_OUT_DIR:
        Full path to directory containing the *.out.bam files.
        Example:
            "/sci/.../nik_data_project/trimmed/STAR_out"
"""

# EDIT THIS to your STAR output directory:
BASE_STAR_OUT_DIR = "FULLPATH/TO/STAR_OUT_DIR"  # should contain *.out.bam
if not BASE_STAR_OUT_DIR.endswith("/"):
    BASE_STAR_OUT_DIR += "/"

os.popen("module load samtools").read()

path = BASE_STAR_OUT_DIR
files = os.popen(f"ls {path}*.out.bam").read().split("\n")[:-1]


def samToBam(sam: str) -> None:
    """Convert SAM → BAM."""
    name = sam.split(".sam")[0]
    bam = os.popen(f"samtools view -b -o {name}.bam {sam}").read()
    print(bam)


def BamToSam(bam: str) -> None:
    """Convert BAM → SAM."""
    name = bam.split(".bam")[0]
    sam = os.popen(f"samtools view -h {bam} -o {name}.sam").read()
    print(sam)


def checkStats(bam: str) -> None:
    """Run samtools flagstat on a BAM."""
    stats = os.popen(f"samtools flagstat {bam}").read()
    print(stats)


def bamToSorted(bam: str) -> None:
    """Sort BAM → <name>.sorted.bam."""
    name = bam.split(".bam")[0]
    sorted_bam = os.popen(f"samtools sort -@ 8 -o {name}.sorted.bam {bam}").read()
    print(sorted_bam)


def indexBam(sorted_bam: str) -> None:
    """Index BAM."""
    index = os.popen(f"samtools index {sorted_bam}").read()
    print(index)


# Example batch loop – currently just BAM → SAM for all files.
for bam in files:
    if not bam:
        continue
    name = bam.split(".bam")[0]

    # Example other steps if you want them:
    # checkStats(bam)
    # bamToSorted(bam)
    # indexBam(f"{name}.sorted.bam")

    BamToSam(bam)

print("end")
