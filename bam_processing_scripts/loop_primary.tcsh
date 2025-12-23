#!/bin/tcsh
###############################################################################
# loop_primary.tcsh
# Author: Noa Gilad
#
# Description:
#   Interactive wrapper to run primary_alignment_filter.py on multiple SAM/BAM
#   files. For each file, submits a SLURM job that separates primary vs
#   non-primary alignments.
#
# Usage:
#   Run inside a directory containing SAM/BAM files (adjust FILE_SUFFIX if
#   needed). The script will prompt for:
#       - A specific sample prefix
#       - Or "ALL" for all matching files
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - Python 3 with primary_alignment_filter.py accessible
#
# Paths requiring user configuration:
#   LOGS          → directory where SLURM log files will be written
#   PYTHON_SCRIPT → full path to primary_alignment_filter.py
###############################################################################

set USAGE = "Usage: loop_orgsep: an interactive script, will ask for the needed parameters"
set tooMany = $1

if ($tooMany != "") then
    echo $USAGE
    exit 1;
endif

### PARAMETERS ###
# Set working directory as the running folder
set SAM_DIR = `pwd`
set LOGS = "/FULLPATH/TO/LOGS"  # Specified SLURM log directory
set NUMT = 16
set MEM = "32G"


# Path to the updated Python script
set PYTHON_SCRIPT = "/FULLPATH/TO/primary_alignment_filter.py"

# Interactive sample selection
set s = ""

while ($s == "")
    echo ""
    echo "*********** Please enter a unique sample prefix for SAM, type ALL for all samples, or type exit to finish ***************"
    set s = "$<"
    if ($s == "exit") then
        continue
    else if ($s == "") then
        continue
    else if ($s == "ALL") then
        # Select all BAM files in the current directory
        set bam_files = `find ${SAM_DIR} -maxdepth 1 -name "*_outputAligned.sortedByCoord.out.bam"`
    else
        # Select BAM files matching a specific prefix
        set bam_files = `find ${SAM_DIR} -maxdepth 1 -name "${s}*_outputAligned.sortedByCoord.out.bam"`
        if (${#bam_files} == 0) then
            echo "No BAM files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through BAM files and submit jobs using sbatch
    foreach bam_file ($bam_files)
        set base = `basename ${bam_file} _Caenorhabditis_elegans_miRNA-Seq_outputAligned.sortedByCoord.out.bam`
        set output_log = "${LOGS}/primarysep_${base}.txt"

        echo ""
        echo "********** Separating BAM by alignment type for sample $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              python3 '${PYTHON_SCRIPT}' '${SAM_DIR}'/'${base}'_celegans.sam ; \
              samtools view -hbS '${SAM_DIR}'/'${base}'_not_primary.sam | samtools sort -@ '${NUMT}' -o '${SAM_DIR}'/'${base}'_not_primary.bam - ; \
              samtools index '${SAM_DIR}'/'${base}'_not_primary.bam ; \
              samtools view -hbS '${SAM_DIR}'/'${base}'_cel.sam | samtools sort -@ '${NUMT}' -o '${SAM_DIR}'/'${base}'_cel.bam - ; \
              samtools index '${SAM_DIR}'/'${base}'_cel.bam ; \
              # Keep C. elegans SAM as SAM without conversion \
              rm -f '${SAM_DIR}'/'${base}'_not_primary.sam '${SAM_DIR}'/'${base}'_cel.sam ;' \
              >! tmp_primary_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${SAM_DIR}/jobid_${base} | cut -d ' ' -f 4`

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=1-0 -o ${output_log} -d afterok:$alid tmp_primary_${base}) >&! ${SAM_DIR}/jobid_primary_${base}


        # Wait briefly to avoid overwhelming the scheduler
        sleep 1
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_primary_*
echo "All BAM files processed and submitted successfully!"
