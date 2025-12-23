#!/bin/tcsh
###############################################################################
# UMI Deduplication (umi_tools)
# Author: Noa Gilad
#
# Description:
#   Interactive script for deduplicating BAM files using umi_tools.
#   The script identifies unique molecules using UMIs embedded in read names,
#   and submits each sample as a SLURM job with sbatch.
#
# Usage:
#   Run inside a directory containing STAR-aligned BAM files named:
#       SAMPLEAligned.sortedByCoord.out.bam
#
#   The script will prompt for:
#       - A specific sample prefix
#       - OR "ALL" to process every BAM file
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - umi_tools installed via conda (umi-tools-env)
#
# Paths requiring user configuration:
#   LOGS     → Path to SLURM log directory
#   JOB_DIR  → Directory containing STAR job ID files (for dependencies)
###############################################################################

set USAGE = "Usage: umi deduplication: an interactive script, will ask for the needed parameters"
set tooMany = $1

if ($tooMany != "") then
    echo $USAGE
    exit 1;
endif

### PARAMETERS ###
# Set working directory as the running folder
set DIR = `pwd`
set LOGS = "/FULLPATH/TO/LOGS"  # Specified SLURM log directory
set NUMT = 16
set MEM = "128G"

echo "*********** Please enter job DIR ***************"
set JOB_DIR = "$<"

# Interactive sample selection
set s = ""

while ($s == "")
    echo ""
    echo "*********** Please enter a unique sample prefix for bam, type ALL for all samples, or type exit to finish ***************"
    set s = "$<"
    if ($s == "exit") then
        continue
    else if ($s == "") then
        continue
    else if ($s == "ALL") then
        # Select all bam files in the current directory
        set bam_files = `find ${DIR} -maxdepth 1 -name "*bam"`
    else
        # Select bam files matching a specific prefix
        set bam_files = `find ${DIR} -maxdepth 1 -name "${s}*.bam"`
        if (${#bam_files} == 0) then
            echo "No bam files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through bam files and submit jobs using sbatch
    foreach bam_file ($bam_files)
        set base = `basename ${bam_file} Aligned.sortedByCoord.out.bam`
        set output_log = "${LOGS}/umi_dedup_${base}.txt"

        echo ""
        echo "********** deduplicating (UMI) $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              module load spack miniconda3 ; \
              conda activate umi-tools-env ; \
              umi_tools dedup -I '${DIR}'/'${base}'Aligned.sortedByCoord.out.bam --paired --output-stats=dedup_stats -S  '${DIR}'/'${base}'_dedup.bam ;' \
              >! tmp_umi_dedup_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${JOB_DIR}/jobid_star_${base} | cut -d ' ' -f 4`

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=2-0 -o ${output_log} -d afterok:$alid tmp_umi_dedup_${base}) >&! ${DIR}/jobid_umi_dedup_${base}


        # Wait briefly to avoid overwhelming the scheduler
        sleep 2
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_umi_dedup_*
echo "All bam files processed and submitted successfully!"
