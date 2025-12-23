#!/usr/bin/env python3
import os
import pandas as pd

"""
rsem_merge_expected_counts.py

Description:
    Merge RSEM per-sample *.genes.results files into a single expected-counts
    matrix and a gene list.

Inputs:
    - RSEM gene-level outputs:
        * <sample>.<org>.genes.results
      All located under:
        BASE_PROJECT_DIR / <date> / "RSEM_output"

Outputs (written into the same RSEM_output directory):
    - genes<date>.txt
        One gene ID (row name) per line.
    - RSEM_matrix_<date>.csv
        Matrix with genes as rows and samples as columns, containing
        the "expected_count" values.

Parameters:
    - DATE: experiment date string used in folder names and output filenames.
    - ORG:  organism short name used in RSEM output filenames.

Paths requiring user configuration:
    - BASE_PROJECT_DIR:
        Full path to the project root that contains <DATE>/RSEM_output.
"""

DATE = "25_1_2022"
ORG = "ce11"

# EDIT THIS to your project root:
BASE_PROJECT_DIR = "FULLPATH/TO/PROJECT_ROOT"  # e.g. /sci/.../project

RSEM_DIR = os.path.join(BASE_PROJECT_DIR, DATE, "RSEM_output")

files_list = (
    os.popen(f'ls {RSEM_DIR}/*.{ORG}.genes.results')
    .read()
    .strip()
    .split("\n")
)

data_frames = []
col_names = []

for file in files_list:
    # Assumes filenames contain ".../RSEM_output/<sample>_R1_001.<org>.genes.results"
    name = file.split("RSEM_output/")[-1].split("_R1_001")[0]
    col_names.append(name)

    df = pd.read_table(file, sep="\t", index_col=0, header=0)
    df = df[["expected_count"]]
    data_frames.append(df)

mat = pd.concat(data_frames, ignore_index=False, axis=1)
mat.columns = col_names

gene_names = mat.index

genes_out = os.path.join(RSEM_DIR, f"genes{DATE}.txt")
matrix_out = os.path.join(RSEM_DIR, f"RSEM_matrix_{DATE}.csv")

with open(genes_out, "w+") as gn:
    for gene in gene_names:
        gn.write(f"{gene}\n")

mat.to_csv(matrix_out, index=False)
