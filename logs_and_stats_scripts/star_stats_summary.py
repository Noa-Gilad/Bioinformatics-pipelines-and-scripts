#!/usr/bin/env python3
###############################################################################
# STAR Alignment Log Summary
# Author: Noa Gilad
#
# Description:
#   Summarizes STAR "Log.final.out" files in a given directory and writes a
#   table with:
#       - total reads
#       - uniquely aligned reads
#       - multi-aligned reads
#       - not aligned reads
#   For aligned / not aligned counts, percentages of total reads are added in
#   the form "num (per%)".
#
# Usage:
#   python star_stats_summary.py <log_dir> <output_path> [filter_str]
#
#   Example:
#     python star_stats_summary.py \
#         /PATH/TO/STAR_out \
#         /PATH/TO/STAR_out \
#         trimmed
#
#   This will create:
#       star_trimmed_summary.txt
#   in <output_path>.
###############################################################################

import os
import sys
import subprocess


def summarize_star_logs(log_dir, output_path, filter_str="trimmed"):
    """
    Summarize STAR alignment logs in 'log_dir' and write the output to
    'output_path'.

    Looks for any file ending with "Log.final.out" and extracts mapping
    statistics.

    For uniquely aligned, multi-aligned, and not aligned values, the script
    adds the percentage of total reads in the format: "num (per%)".

    Args:
        log_dir (str): Directory containing STAR log files.
        output_path (str): Directory to save the summary file.
        filter_str (str, optional): Label (e.g. 'trimmed').
                                    Default is 'trimmed'.
    """
    log_files = [f for f in os.listdir(log_dir) if f.endswith("Log.final.out")]
    if not log_files:
        print(f"No STAR log files ending with 'Log.final.out' found in {log_dir}")
        return

    summary_data = []

    for f in log_files:
        # Sample name is everything before "Log.final.out"
        sample_name = f.replace("Log.final.out", "")
        input_log = os.path.join(log_dir, f)

        # Extract total reads
        cmd_total = (
            f"grep 'Number of input reads' '{input_log}' | "
            "awk -F'|' '{print $2}'"
        )
        total_reads = subprocess.run(
            cmd_total,
            shell=True,
            capture_output=True,
            text=True
        ).stdout.strip()

        # Extract uniquely aligned reads
        cmd_unique = (
            f"grep 'Uniquely mapped reads number' '{input_log}' | "
            "awk -F'|' '{print $2}'"
        )
        uniquely_aligned = subprocess.run(
            cmd_unique,
            shell=True,
            capture_output=True,
            text=True
        ).stdout.strip()

        # Extract multi-aligned reads
        cmd_multi = (
            f"grep 'Number of reads mapped to multiple loci' '{input_log}' | "
            "awk -F'|' '{print $2}'"
        )
        multi_aligned = subprocess.run(
            cmd_multi,
            shell=True,
            capture_output=True,
            text=True
        ).stdout.strip()

        # Convert values to integers and compute not aligned
        if total_reads and uniquely_aligned and multi_aligned:
            total_reads = int(total_reads)
            uniquely_aligned = int(uniquely_aligned)
            multi_aligned = int(multi_aligned)
            not_aligned = total_reads - uniquely_aligned - multi_aligned
        else:
            print(f"Warning: Missing values in {input_log}. Skipping.")
            continue

        summary_data.append(
            (sample_name, total_reads,
             uniquely_aligned, multi_aligned, not_aligned)
        )

    # Sort the summary data by sample name
    summary_data.sort(key=lambda x: x[0])

    summary_filename = f"star_{filter_str}_summary.txt"
    summary_path = os.path.join(output_path, summary_filename)

    with open(summary_path, 'w') as summary_file:
        summary_file.write(f"#Summary of STAR ({filter_str})\n")
        summary_file.write(
            "#sample\tTotal_reads\tUniquely_aligned\t"
            "Multi_aligned\tNot_aligned\n"
        )
        for row in summary_data:
            sample_name, total, unique, multi, not_aligned = row
            unique_str = f"{unique} ({unique / total * 100:.1f}%)"
            multi_str = f"{multi} ({multi / total * 100:.1f}%)"
            not_aligned_str = f"{not_aligned} ({not_aligned / total * 100:.1f}%)"
            summary_file.write(
                f"{sample_name}\t{total}\t"
                f"{unique_str}\t{multi_str}\t{not_aligned_str}\n"
            )

    print(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    """
    Usage:
      python star_stats_summary.py <log_dir> <output_path> [filter_str]

    Example:
      python star_stats_summary.py /PATH/TO/STAR_out /PATH/TO/STAR_out trimmed
    """
    if len(sys.argv) < 3:
        print(
            "Usage: python star_stats_summary.py <log_dir> <output_path> [filter_str]"
        )
        sys.exit(1)

    log_dir = sys.argv[1]
    output_path = sys.argv[2]
    filter_str = sys.argv[3] if len(sys.argv) > 3 else "trimmed"

    summarize_star_logs(log_dir, output_path, filter_str)
