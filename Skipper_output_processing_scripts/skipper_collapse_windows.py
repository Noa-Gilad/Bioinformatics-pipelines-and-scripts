#!/usr/bin/env python3
###############################################################################
# Skipper per-gene collapse of window-level output
# Author: Noa Gilad
#
# Description:
#   Takes Skipper window-level output (multiple rows per gene_name) and
#   collapses it to a single row per gene:
#       - windows_coor: all window coordinates for that gene
#       - peaks_coor: merged intervals (non-overlapping) spanning windows
#       - selected annotation/statistical columns collapsed per gene
#
# Usage:
#   python skipper_collapse_windows.py <input.tsv[.gz]> <output.tsv>
#
# Input:
#   - Skipper windows table with at least:
#       gene_name, start, end, chr/chr-like, strand, stats columns...
#
# Output:
#   - TSV with one row per gene_name, containing window/peak coordinates and
#     aggregated columns.
#
# Notes:
#   - DROP_COLS, SINGLE, and ORDER control which columns are kept and how.
###############################################################################

import sys
import pandas as pd

# --- config ---
SINGLE = {
    "gene_name", "chr", "strand",
    "gene_id", "transcript_ids", "gene_type_top"   # <- newly forced single
}
SEP = ", "

DROP_COLS = {
    "transcript_type_top", "gene_types", "transcript_types",
    "feature_type_top", "feature_types", "feature_id", "feature_bin",
    "score", "chrom"
}

ORDER = [
    "gene_name", "chr", "strand",
    "peaks_coor", "windows_coor",
    "name", "score", "gc", "gc_bin", "chrom",
    "gene_id", "transcript_ids", "gene_type_top",
    "input_sum", "clip_sum", "enrichment_n",
    "enrichment_l2or_min", "enrichment_l2or_mean", "enrichment_l2or_max",
    "p_max", "p_min", "q_max", "q_min"
]

def merge_intervals(starts, ends):
    ivals = sorted(zip(starts, ends))
    merged = []
    for s, e in ivals:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)
    return SEP.join(f"{s}:{e}" for s, e in merged)

def collapse(group: pd.DataFrame) -> pd.Series:
    g = group.sort_values(["start", "end"])
    out = {}
    out["windows_coor"] = SEP.join((g["start"].astype(str)+":"+g["end"].astype(str)).tolist())
    out["peaks_coor"]   = merge_intervals(g["start"].tolist(), g["end"].tolist())

    for col in g.columns:
        if col in {"start", "end", "windows_coor", "peaks_coor"}:
            continue
        vals = g[col].astype(str).tolist()
        out[col] = vals[0] if col in SINGLE else SEP.join(vals)
    return pd.Series(out)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <input.tsv[.gz]> <output.tsv>\n")
        sys.exit(1)

    IN, OUT = sys.argv[1], sys.argv[2]
    df = pd.read_csv(IN, sep="\t", compression="infer")

    df = df.drop(columns=[c for c in DROP_COLS if c in df.columns], errors="ignore")
    df["start"] = pd.to_numeric(df["start"])
    df["end"]   = pd.to_numeric(df["end"])

    res = (df.groupby("gene_name", as_index=False)
             .apply(collapse)
             .reset_index(drop=True))

    cols = [c for c in ORDER if c in res.columns]
    res = res[cols + [c for c in res.columns if c not in cols]]

    res.to_csv(OUT, sep="\t", index=False)
