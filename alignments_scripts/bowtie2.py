#!/usr/bin/env python3
"""
bowtie2_batch_align.py

Submit a Bowtie2 alignment job per FASTQ file using SLURM.

For each *.fastq in a given experiment folder, this script:
  - Builds an output SAM path under bowtie2_out/
  - Aligns reads with a given Bowtie2 index (mi21U)
  - Submits the command via `sbatch`

Usage
-----
    python bowtie2_batch_align.py

Configuration
-------------
Edit the values in the CONFIG section below:

- BASE_PROJECT_DIR : base directory containing the experiment subfolder
- EXPERIMENT       : experiment subfolder name (e.g. "trimmed/")
- BOWTIE2_INDEX_DIR: directory containing the Bowtie2 index prefix
- LOGS_DIR         : directory for SLURM log files

Notes
-----
- This script assumes a single-end layout (parameter -U).
- The exact index prefix (here "mi21U") can be changed if needed.
"""

import os

# --------------------- CONFIG --------------------- #
BASE_PROJECT_DIR = "/FULLPATH/TO/PROJECT"
EXPERIMENT = "SUBDIRECTORY/"
BOWTIE2_INDEX_DIR = "/FULLPATH/TO/INDEX"
LOGS_DIR = "/FULLPATH/TO/LOGS"
INDEX_NAME = "mi21U"
# -------------------------------------------------- #

os.popen("module load bowtie2").read()

experiment = EXPERIMENT
path = f"{BASE_PROJECT_DIR}/{experiment}"
files = os.popen(f'ls {path}*.fastq').read().split("\n")[:-1]

for file in files:
    if not file:
        continue
    print(file)
    name = file.split(experiment)[-1].split(".fastq")[0]
    out_file = f"{path}/bowtie2_out/{name}"
    ref = f"{BOWTIE2_INDEX_DIR}/{INDEX_NAME}"
    command = f"bowtie2 -x {ref} -U {file} -S {out_file}.sam"

    sbatch_cmd = (
        f'sbatch --mem=128g --time=2-0 '
        f'--job-name=align_{name} '
        f'-o {LOGS_DIR}/bowtie2.txt '
        f'--wrap="{command}"'
    )
    out = os.popen(sbatch_cmd).read()
    print(out)

print("end")
