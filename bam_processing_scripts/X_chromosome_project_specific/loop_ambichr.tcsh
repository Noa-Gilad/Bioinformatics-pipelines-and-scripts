#!/bin/tcsh
###############################################################################
# loop_ambichr.tcsh
# Author: Noa Gilad
#
# Description:
#   Interactive wrapper to run ambi_chr_sep.py on multiple C. elegans SAM files.
#   For each SAM, submits a SLURM job that splits reads into:
#     - X-only alignments
#     - X + autosome (ambiguous chromosome) alignments.
#
# Usage:
#   Run inside a directory containing the relevant SAM/BAM files and follow
#   the prompts:
#       - sample prefix
#       - or "ALL" to process all matching files.
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - Python 3 with ambi_chr_sep.py available
#
# Paths requiring user configuration:
#   LOGS          → directory where SLURM log files will be written
#   PYTHON_SCRIPT → full path to ambi_chr_sep.py
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
set LOGS = "/FULLPATH/TO/LOGS"               # Specified SLURM log directory
set PYTHON_SCRIPT = "/FULLPATH/TO/ambi_chr_sep.py"
set NUMT = 16
set MEM = "32G"

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
        set output_log = "${LOGS}/chrsep_${base}.txt"

        echo ""
        echo "********** Separating BAM by chr for sample $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              python3 '${PYTHON_SCRIPT}' '${SAM_DIR}'/'${base}'_celegans.sam ; \
              samtools view -hbS '${SAM_DIR}'/'${base}'_cel_ambiguous_chr.sam | samtools sort -@ '${NUMT}' -o '${SAM_DIR}'/'${base}'_cel_ambiguous_chr.bam - ; \
              samtools index '${SAM_DIR}'/'${base}'_cel_ambiguous_chr.bam ; \
              # Keep C. elegans SAM as SAM without conversion \
              rm -f '${SAM_DIR}'/'${base}'_cel_ambiguous_chr.sam ;' \
              >! tmp_ambi_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${SAM_DIR}/jobid_${base} | cut -d ' ' -f 4`

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=1-0 -o ${output_log} -d afterok:$alid tmp_ambi_${base}) >&! ${SAM_DIR}/jobid_cel_chr_${base}


        # Wait briefly to avoid overwhelming the scheduler
        sleep 1
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_ambi_*
echo "All BAM files processed and submitted successfully!"
