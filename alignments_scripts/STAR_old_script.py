#!/usr/bin/env python3
"""
STAR_old_script.py
Author: Noa Gilad

Description:
  Submit STAR alignment jobs (one per FASTQ) via SLURM.

Before running (user config):
  - Set experiment to your sub folder name (example shown)
  - Set BASE_DIR to your project directory
  - Set STAR_INDEX to your STAR genomeDir
  - Set LOG_FILE to where you want SLURM logs to go
"""

import os

os.popen('module load STAR').read()

experiment = "trimmed/"  # example
BASE_DIR = "/FULLPATH/TO/PROJECT_DIR/"  # e.g. "/data/my_project/"
path = f"{BASE_DIR}{experiment}"

files = os.popen(f'ls {path}*.fastq').read().split("\n")[:-1]

STAR_INDEX = "/FULLPATH/TO/STAR_INDEX/"  # genomeDir directory
LOG_FILE = "/FULLPATH/TO/LOGS/STAR.txt"

for file in files:
    print(file)
    name = file.split(experiment)[-1].split(".fastq")[0]
    out_file = f"{path}/STAR_out/{name}"

    ref = STAR_INDEX
    command = (
        f"STAR --runMode alignReads --genomeDir {ref} "
        f"--outSAMtype BAM SortedByCoordinate --readFilesIn {file} "
        f"--runThreadN 12 --outFileNamePrefix {out_file} "
        f"--outFilterMultimapNmax 1000 --outSAMprimaryFlag AllBestScore"
    )

    out = os.popen(
        f'sbatch --mem=64g --time=1-0 --job-name=STAR -o {LOG_FILE} --wrap="{command}"'
    ).read()
    print(out)

print("end")
