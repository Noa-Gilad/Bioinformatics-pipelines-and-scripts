#!/usr/bin/env python3
"""
RSEM.py
Author: Noa Gilad

Description:
  Submit RSEM expression quantification jobs (one per FASTQ) via SLURM.

Before running (user config):
  - Set experiment to your subfolder name (example shown)
  - Set BASE_DIR to your project directory (NOT lab-specific)
  - Set org to your organism build (example: "ce11")
  - Set RSEM_REF to your RSEM reference prefix
  - Set LOG_FILE to where you want SLURM logs to go
"""

import os

os.popen('module load rsem')

experiment = "25_1_2022/"  # example
BASE_DIR = "/FULLPATH/TO/PROJECT_DIR/"  # e.g. "/data/my_project/"
path = f"{BASE_DIR}{experiment}"

files = os.popen(f'ls {path}*.fastq').read().split("\n")[:-1]

org = "ce11"  # organism build, e.g. "ce11"
RSEM_REF = f"/FULLPATH/TO/RefSeq_RSEM/{org}/RefSeq_RSEM"  # adjust to your needs
LOG_FILE = "/FULLPATH/TO/LOGS/RSEM.txt"

for file in files:
    print(file)
    name = file.split(experiment)[-1].split(".fastq")[0]
    out_file = f"{path}/RSEM_output/{name}.{org}"

    command = f"rsem-calculate-expression --num-threads 1 --bowtie2 --estimate-rspd {file} {RSEM_REF} {out_file}"
    out = os.popen(
        f'sbatch --mem=3g --time=2-0 --job-name=RSEM -o {LOG_FILE} --wrap="{command}"'
    ).read()
    print(out)

print("end")
