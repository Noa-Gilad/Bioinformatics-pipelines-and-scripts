#!/usr/bin/env python3
"""
x-chromosome stats plotting

Parse per-sample stats text files for the X-chromosome project and generate
summary plots for mapping percentages, read counts, and small RNA biotypes.

Expected input:
    <BASE_STATS_DIR>/<EXPERIMENT_SUBDIR>/<sample>_stats.txt

Outputs multiple PNG figures into an output folder.
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# ---- configuration ----
# Base directory containing the stats files
BASE_STATS_DIR = "/FULLPATH/TO/x_chromosome_project/nik_data_project"
# Subdirectory for a specific run (relative to BASE_STATS_DIR)
EXPERIMENT_SUBDIR = "trimmed/STAR_out/"

# Compose full path
path = f"{BASE_STATS_DIR}/{EXPERIMENT_SUBDIR}"

# List all *_stats.txt files
files = os.popen(f'ls {path}*_stats.txt').read().split("\n")[:-1]

# Storage for statistics
samples = []
mapped_to_celegans_pct = []
mapped_to_celegans_reads = []
biotype_reads = {}
biotype_genes = {}

# Storage for total small RNA reads per sample
small_rna_reads = []
small_rna_genes = []

# Small RNA biotypes
small_rna_types = ["snRNA", "snoRNA", "tRNA", "miRNA", "piRNA"]

# Loop over each file
for file in files:
    sample_name = file.split(EXPERIMENT_SUBDIR)[-1].replace("_stats.txt", "")
    samples.append(sample_name)

    with open(file, "r") as f:
        content = f.read()

    # Extract percentage of reads mapped to c.elegans
    pct_match = re.search(r'Percentage of reads mapped to c\.elegans only and has the highest score from all mapped reads\t([\d.]+)', content)
    mapped_to_celegans_pct.append(float(pct_match.group(1)) if pct_match else 0)

    # Extract number of reads mapped to c.elegans
    reads_match = re.search(r'Mapped reads to c\.elegans only and has the highest score\t(\d+)', content)
    mapped_to_celegans_reads.append(int(reads_match.group(1)) if reads_match else 0)

    # Track small RNA totals for this sample
    small_rna_count = 0
    small_rna_gene_count = 0

    # Parse small RNA reads and genes
    for biotype in small_rna_types:
        # Reads per biotype
        reads_match = re.search(f'{biotype}, (\\d+)', content)
        biotype_count = int(reads_match.group(1)) if reads_match else 0
        biotype_reads.setdefault(biotype, []).append(biotype_count)

        # Add to small RNA total counter
        small_rna_count += biotype_count

        # Genes with more than 2 reads per biotype
        genes_match = re.search(f'Genes With More Than 2 Reads Per Biotype\n.*{biotype}, (\\d+)', content, re.DOTALL)
        gene_count = int(genes_match.group(1)) if genes_match else 0
        biotype_genes.setdefault(biotype, []).append(gene_count)

        # Add to small RNA gene counter
        small_rna_gene_count += gene_count

    # Store total small RNA counts for this sample
    small_rna_reads.append(small_rna_count)
    small_rna_genes.append(small_rna_gene_count)


# Function to calculate average excluding outliers
def calculate_average_without_outliers(data):
    z_scores = np.abs(stats.zscore(data))
    filtered_data = [x for x, z in zip(data, z_scores) if z < 3]
    return np.mean(filtered_data)


# Function to plot histograms with value labels
def plot_histogram(x, y, title, xlabel, ylabel, filename, show_average=True):
    plt.figure(figsize=(15, 5))
    bars = plt.bar(x, y, color="skyblue")

    # Calculate and add average (excluding outliers)
    if show_average:
        avg = calculate_average_without_outliers(y)
        title += f'\nAverage (excluding outliers): {avg:.2f}'

    plt.title(title)
    plt.xticks(rotation=90, fontsize=8)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2., height,
                 f'{height:.1f}',
                 ha='center', va='bottom', fontsize=8, rotation=90)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


# Function to plot a distribution histogram for reads
def plot_reads_distribution(mapped_to_celegans_reads):
    plt.figure(figsize=(10, 5))

    # Create bins from 0 to max value with 5000 interval
    max_reads = max(mapped_to_celegans_reads)
    bins = list(range(0, int(max_reads) + 10001, 10000))

    plt.hist(mapped_to_celegans_reads, bins=bins, edgecolor='black', color="lightcoral")
    plt.title('Distribution of Reads per Sample')
    plt.xlabel('Number of Reads')
    plt.ylabel('Number of Samples')
    plt.xticks(bins, rotation=90)
    plt.tight_layout()
    plt.savefig('biotype_analysis_plots/reads_distribution.png')
    plt.close()


# Create output directory
os.makedirs('biotype_analysis_plots', exist_ok=True)

# C. elegans mapping plots
plot_histogram(samples, mapped_to_celegans_pct,
               "Percentage of Reads Mapped to C. elegans Only (all alignments with the highest score)",
               "Sample", "Percentage (%)",
               "biotype_analysis_plots/mapped_to_celegans_pct.png",
               show_average=False)  # Don't calculate average for percentages

plot_histogram(samples, mapped_to_celegans_reads,
               "Mapped Reads to C. elegans (all alignments with the highest score)",
               "Sample", "Number of Reads",
               "biotype_analysis_plots/mapped_to_celegans_reads.png")

# Small RNA-specific plots
for biotype in small_rna_types:
    # Reads per small RNA biotype
    plot_histogram(samples, biotype_reads[biotype],
                   f"Reads Mapped to {biotype}",
                   "Sample", "Number of Reads",
                   f"biotype_analysis_plots/{biotype}_reads.png")

    # Genes with more than 2 reads per small RNA biotype
    plot_histogram(samples, biotype_genes[biotype],
                   f"Genes with >2 Reads in {biotype}",
                   "Sample", "Number of Genes",
                   f"biotype_analysis_plots/{biotype}_genes.png")

# New histogram for total reads distribution
plot_reads_distribution(mapped_to_celegans_reads)

# Small RNA reads histogram
plot_histogram(samples, small_rna_reads,
               "Total Small RNA Reads per Sample",
               "Sample", "Number of Reads",
               "biotype_analysis_plots/small_rna_reads.png")

# Small RNA genes histogram
plot_histogram(samples, small_rna_genes,
               "Genes with >2 Reads in Small RNAs",
               "Sample", "Number of Genes",
               "biotype_analysis_plots/small_rna_genes.png")

print("✅ All graphs saved successfully!")
