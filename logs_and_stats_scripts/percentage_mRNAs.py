#!/usr/bin/env python3
###############################################################################
# Percentage of mRNA Reads in Count Table
# Author: Noa Gilad
#
# Description:
#   Given:
#       1) A TSV file with a list of mRNA gene IDs (column "WB" - match Wormbase datset)
#       2) A featureCounts/HTSeq count table (TSV) with a "Geneid" column
#   this script calculates what percentage of total counts come from the
#   mRNA gene list.
#
# Usage:
#   1. Set MRNA_LIST_TSV and COUNT_TABLE_TSV in the __main__ section.
#   2. Run:
#          python percentage_mRNAs.py
#   3. The script prints the percentage of reads assigned to mRNA transcripts.
###############################################################################

import pandas as pd


def read_gene_list(file_path):
    """
    Read a TSV file containing a "WB" column and return a set of gene IDs.
    """
    df = pd.read_csv(file_path, sep='\t')
    return set(df['WB'])


def check_in_list(genes, file_path):
    """
    Calculate the percentage of counts assigned to a given gene set.

    Args:
        genes (set): A set of gene IDs (strings).
        file_path (str): Path to a count table TSV with a "Geneid" column.

    Returns:
        float: Percentage of total counts that belong to 'genes'.
    """
    df = pd.read_csv(file_path, sep='\t', skiprows=1)

    total = 0
    mrna_counts = 0

    for _, row in df.iterrows():
        # Sum over all samples (columns 7 onward)
        sum_row = row.iloc[6:].sum()
        total += sum_row

        if row['Geneid'] in genes:
            mrna_counts += sum_row

    print(f"total counts (all genes): {total}")
    percentage = mrna_counts / total * 100
    return percentage


if __name__ == '__main__':
    # Paths to configure:
    #   MRNA_LIST_TSV  – TSV with a "WB" column listing mRNA genes
    #   COUNT_TABLE_TSV – count table from featureCounts / HTSeq
    MRNA_LIST_TSV = "PATH/TO/mRNAs.tsv"
    COUNT_TABLE_TSV = "PATH/TO/count.out"

    gene_list = read_gene_list(MRNA_LIST_TSV)
    percent = check_in_list(gene_list, COUNT_TABLE_TSV)
    print(f"the percentage is {percent}")
