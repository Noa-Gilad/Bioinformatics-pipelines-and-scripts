#!/usr/bin/env python3
###############################################################################
# Cutadapt Log Parser and Poly-Adapter Statistics
# Author: Noa Gilad
#
# Description:
#   Parses Cutadapt log files and extracts:
#       - Percentage of reads with adapters
#       - Percentages of polyA / polyT / polyG / polyC trimming
#   For each log file, creates a bar plot and writes a summary table.
#
# Usage:
#   1. Place this script in a directory with Cutadapt log files named like:
#          cutadapt_no_trimming_<SAMPLE>.txt
#   2. Run the script (e.g. from that directory):
#          python cutadapt_log_parser.py
#   3. The script will:
#          - Create per-sample PNG plots
#          - Create a summary file: cutadapt_summary.txt
#
# Notes:
#   - You can change the output directory in the __main__ section below.
###############################################################################

import os
import re
import matplotlib.pyplot as plt


def parse_cutadapt_log(log_content):
    """
    Parse a Cutadapt log file and extract percentages.

    The log file is expected to contain two parts:
      1. The original summary (with "Total reads processed" and
         "Reads with adapters").
      2. Poly adapter sections (polyA, polyT, polyG, polyC).

    For the poly sections, the trimmed counts are extracted using a non-greedy
    match so that the first "Trimmed:" after the header is used.

    Args:
        log_content (str): Content of the Cutadapt log file.

    Returns:
        dict: Percentages of adapters and poly types relative to total reads.
              Keys: 'adapter', 'polyA', 'polyT', 'polyG', 'polyC'.
    """
    # Extract total reads from the first summary section
    total_reads_match = re.search(
        r'Total reads processed:\s*([0-9,]+)',
        log_content
    )
    if total_reads_match:
        total_reads = int(total_reads_match.group(1).replace(',', ''))
    else:
        raise ValueError("Could not find total reads in the log file")

    # Extract adapter trimming percentage from the first summary section
    adapter_match = re.search(
        r'Reads with adapters:\s*[0-9,]+\s+\(([0-9.]+)%\)',
        log_content
    )
    adapter_percent = float(adapter_match.group(1)) if adapter_match else 0.0

    # Initialize dictionary with adapter percentage and zero for each poly type
    poly_percentages = {
        'adapter': adapter_percent,
        'polyA': 0.0,
        'polyT': 0.0,
        'polyG': 0.0,
        'polyC': 0.0,
    }

    # For each poly type, use a non-greedy match (.*?) so that we capture
    # the first occurrence of "Trimmed:" in that section.
    poly_types = ['polyA', 'polyT', 'polyG', 'polyC']
    for poly_type in poly_types:
        pattern = (
            rf'===\s*Adapter\s+{poly_type}\s*===\s*.*?'
            r'Trimmed:\s*([0-9,]+)\s+times'
        )
        poly_match = re.search(pattern, log_content, re.DOTALL)
        if poly_match:
            poly_reads = int(poly_match.group(1).replace(',', ''))
            poly_percentages[poly_type] = (poly_reads / total_reads) * 100.0
        else:
            # If the section is not found, keep it at 0.0%
            poly_percentages[poly_type] = 0.0

    return poly_percentages


def create_bar_plot(analysis_results, output_filename):
    """
    Create a bar plot of adapter and poly-nucleotide trimming percentages.

    Args:
        analysis_results (dict): Dictionary with keys:
            adapter, polyA, polyT, polyG, polyC.
        output_filename (str): Path (with file name) to save the output plot.
    """
    categories = list(analysis_results.keys())
    percentages = list(analysis_results.values())

    plt.figure(figsize=(10, 6))

    # Simple fixed color list; adjust as you like
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    plt.bar(categories, percentages, color=colors)

    plt.title('Adapter and Poly-Nucleotide Trimming Percentages', fontsize=15)
    plt.ylabel('Percentage of Total Reads', fontsize=12)
    plt.ylim(0, 100)
    plt.xticks(rotation=45, ha='right')

    # Add labels above each bar
    for i, v in enumerate(percentages):
        plt.text(i, v + 1, f'{v:.1f}%', ha='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()


def process_log_files(output_dir='.'):
    """
    Process all Cutadapt log files in the current directory that start with
    "cutadapt_no_trimming_".

    For each file:
      - The sample name is extracted (removing "cutadapt_no_trimming_"
        and ".txt").
      - A bar plot is saved as "<sample_name>_stats.png".
      - A summary file "cutadapt_summary.txt" is created in 'output_dir'.

    Args:
        output_dir (str): Directory to save output files.
                          Defaults to current directory.
    """
    os.makedirs(output_dir, exist_ok=True)

    log_files = [
        f for f in os.listdir('.')
        if f.startswith('cutadapt_no_trimming_') and f.endswith('.txt')
    ]

    summary_file_path = os.path.join(output_dir, 'cutadapt_summary.txt')
    with open(summary_file_path, 'w') as summary_file:
        summary_file.write(
            "Sample\tAdapter %\tPolyA %\tPolyT %\tPolyG %\tPolyC %\n"
        )

        for log_file in log_files:
            sample_name = (
                log_file
                .replace('cutadapt_no_trimming_', '')
                .replace('.txt', '')
            )

            try:
                with open(log_file, 'r') as f:
                    log_content = f.read()

                analysis_results = parse_cutadapt_log(log_content)

                output_filename = os.path.join(
                    output_dir,
                    f'{sample_name}_stats.png'
                )
                create_bar_plot(analysis_results, output_filename)

                summary_file.write(
                    f"{sample_name}\t" +
                    "\t".join(
                        f"{analysis_results[poly]:.1f}"
                        for poly in ['adapter', 'polyA', 'polyT', 'polyG', 'polyC']
                    ) +
                    "\n"
                )

                print(f"\nResults for {sample_name}:")
                for poly_type, percentage in analysis_results.items():
                    print(f"{poly_type} Reads: {percentage:.1f}%")

            except Exception as e:
                print(f"Error processing {log_file}: {e}")

    print(f"\nAnalysis complete. Results and plots saved in {output_dir}")


if __name__ == "__main__":
    # Default: write plots and summary into a local subdirectory.
    # You can change this to any folder, e.g. "/PATH/TO/cutadapt/statistics".
    process_log_files("./cutadapt_stats_output")
