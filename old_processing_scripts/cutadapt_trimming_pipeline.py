#!/usr/bin/env python3
import os

"""
cutadapt_trimming_pipeline.py

Description:
    Legacy trimming pipeline using cutadapt for a series of FASTQ files.
    For each input FASTQ, it runs multiple trimming steps (adapter, A-tail,
    spacers, overhang) and writes intermediate and final outputs to an
    output directory.

Current behavior (as written):
    - Loads cutadapt via "module load cutadapt".
    - Lists all *.fastq files under INPUT_DIR.
    - For each file, runs:
        trim_adapter  → trim_A → trim_overhang
      (other helper functions exist but are not used in the final loop).

Usage:
    1. Edit INPUT_DIR and OUTPUT_DIR below.
    2. Make sure cutadapt is installed or loadable via module.
    3. Run the script; it will process all *.fastq files in INPUT_DIR.

Requirements:
    - cutadapt installed and available in PATH or via "module load cutadapt".

Paths requiring user configuration:
    - INPUT_DIR:
        Directory containing the raw *.fastq files.
    - OUTPUT_DIR:
        Directory where trimmed FASTQs should be written.
"""

# EDIT THESE to your environment:
INPUT_DIR = "FULLPATH/TO/RAW_SLICES_DIR/"      # e.g. /sci/.../slices/
OUTPUT_DIR = "FULLPATH/TO/TRIMMED_OUTPUT_DIR/"  # e.g. /sci/.../trimmed/

if not INPUT_DIR.endswith("/"):
    INPUT_DIR += "/"
if not OUTPUT_DIR.endswith("/"):
    OUTPUT_DIR += "/"

os.popen("module load cutadapt").read()

path = INPUT_DIR
files = os.popen(f"ls {path}*.fastq").read().split("\n")[:-1]

output_path = OUTPUT_DIR


def trim_adapter(file, adapter):
    name = file.split(path)[1]
    name = name.split(".fastq")[0]
    os.popen(
        f"cutadapt -a {adapter} --minimum-length=15  --maximum-length=50 "
        f"--quality-cutoff=30 -o {output_path}{name}_tmp1.fastq {file}"
    ).read()


def trim_spacer_1(file, adapter):
    name = file.split(path)[1]
    name = name.split(".fastq")[0]
    os.popen(
        f"cutadapt -a {adapter} --minimum-length=15 "
        f"-o {output_path}{name}_tmp2.fastq {output_path}{name}_tmp1.fastq"
    ).read()


def trim_spacer_2(file, adapter):
    name = file.split(path)[1]
    name = name.split(".fastq")[0]
    os.popen(
        f"cutadapt -a {adapter} --minimum-length=15 "
        f"-o {output_path}{name}_tmp2.fastq {output_path}{name}_tmp1.fastq"
    ).read()


def trim_A(file):
    name = file.split(path)[1]
    name = name.split(".fastq")[0]
    os.popen(
        "cutadapt -a AAAAAAAAAA  --minimum-length=15  "
        f"--maximum-length=50 --quality-cutoff=30 "
        f"-o {output_path}{name}_tmp2.fastq {output_path}{name}_tmp1.fastq"
    ).read()


def trim_overhang(file):
    name = file.split(path)[1]
    name = name.split(".fastq")[0]
    os.popen(
        "cutadapt -u 3  --minimum-length=15  --maximum-length=50 "
        f"--quality-cutoff=30 -o {output_path}{name}_output.fastq "
        f"{output_path}{name}_tmp2.fastq"
    ).read()


def trim_bases(file):
    name = file.split(path)[1]
    name = name.split(".fastq")[0]
    os.popen(
        f"cutadapt -u 3 --minimum-length=15 --maximum-length=50 "
        f"--quality-cutoff=30 -o {output_path}{name}_output.fastq {file}"
    ).read()


# Optional cleanup of intermediates (currently commented out):
# def remove_intermediate_files(file):
#     name = file.split(path)[1]
#     name = name.split(".fastq")[0]
#     os.popen(
#         f"rm {output_path}{name}_tmp1.fastq {output_path}{name}_tmp2.fastq"
#     ).read()


for file in files:
    if not file:
        continue
    trim_adapter(file, "GATCGGAAGAGCACACGTCTGAACTCCAGTCA")
    trim_A(file)
    trim_overhang(file)
    print(f"{file} got trimmed")

print("done")
