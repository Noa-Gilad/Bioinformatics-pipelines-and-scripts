#!/usr/bin/env python3
###############################################################################
# Compute Per-LincRNA Enrichment Scores
# Author: Noa Gilad
#
# Description:
#   Reads a per-window Skipper TSV (gzipped or not) with columns including:
#       gene_name, clip_sum, input_sum, enrichment_l2or_mean
#
#   For each gene_name (lincRNA) it computes:
#       - log2_global_ratio:
#           log2( sum(clip_sum) / sum(input_sum) )
#       - max_window_log2_enrichment:
#           max(enrichment_l2or_mean) across all windows of that gene
#
#   Outputs a gene-level TSV sorted by max_window_log2_enrichment descending.
#
# Usage:
#   python compute_per_linc_scores.py <input.tsv[.gz]> <output.tsv>
#
# Input:
#   <input.tsv[.gz]>  Per-window Skipper table.
#
# Output:
#   <output.tsv>      Gene-level table with:
#                        gene_name
#                        log2_global_ratio
#                        max_window_log2_enrichment
#
# Requirements:
#   - Python 3
#   - pandas
#   - numpy
###############################################################################

import sys
import pandas as pd
import numpy as np

"""
compute_per_linc_scores.py

Reads a per-window TSV (gzipped or not) with columns including:
  gene_name, clip_sum, input_sum, enrichment_l2or_mean

Computes, for each gene_name (lincRNA):
  - log2_global_ratio: log2(sum(clip_sum) / sum(input_sum))
  - max_window_log2_enrichment: max(enrichment_l2or_mean) across windows

Outputs a two-column TSV sorted by max_window_log2_enrichment descending.
"""


def compute_per_linc_enrichment_scores(df: pd.DataFrame) -> pd.DataFrame:
    # 1) Global log2 ratio: log2(sum clip_sum / sum input_sum)
    agg = df.groupby('gene_name').agg({'clip_sum': 'sum', 'input_sum': 'sum'})
    agg['log2_global_ratio'] = np.log2(agg['clip_sum'] / agg['input_sum'])
    # 2) Best window log2 enrichment
    best = (df.groupby('gene_name')['enrichment_l2or_mean']
            .max()
            .rename('max_window_log2_enrichment'))
    # Combine and sort
    res = agg.join(best).reset_index()
    res = res[['gene_name', 'log2_global_ratio', 'max_window_log2_enrichment']]
    res = res.sort_values('max_window_log2_enrichment', ascending=False)
    return res


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <input.tsv.gz> <output.tsv>\n")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_tsv = sys.argv[2]

    # Read the per-window table
    df = pd.read_csv(input_tsv, sep="\t", compression="infer")
    # Ensure numeric
    df['clip_sum'] = pd.to_numeric(df['clip_sum'])
    df['input_sum'] = pd.to_numeric(df['input_sum'])
    df['enrichment_l2or_mean'] = pd.to_numeric(df['enrichment_l2or_mean'])

    # Compute per-linc scores
    result = compute_per_linc_enrichment_scores(df)

    # Write to output
    result.to_csv(output_tsv, sep="\t", index=False)
