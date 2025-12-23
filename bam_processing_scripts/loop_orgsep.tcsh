#!/bin/tcsh
###############################################################################
# loop_orgsep.tcsh
# Author: Noa Gilad
#
# Description:
#   Interactive wrapper to run org_sep.py on multiple BAM files.
#   For each BAM, submits a SLURM job that calls a Python script which splits
#   reads into E. coli, C. elegans, and ambiguous groups.
#
# Usage:
#   Run inside a directory containing BAM files, e.g.:
#       SAMPLE_outputAligned.sortedByCoord.out.bam
#   The script will prompt for:
#       - A specific sample prefix
#       - Or "ALL" for all matching BAMs
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - Python 3 with org_sep.py accessible
#
# Paths requiring user configuration:
#   LOGS          → directory where SLURM log files will be written
#   PYTHON_SCRIPT → full path to org_sep.py
###############################################################################

set USAGE = "Usage: loop_orgsep: an interactive script, will ask for the needed parameters"
set tooMany = $1

if ($tooMany != "") then
    echo $USAGE
    exit 1;
endif

### PARAMETERS ###
# Set working directory as the running folder
set BAM_DIR = `pwd`
set LOGS = "/FULLPATH/TO/LOGS"  # Specified SLURM log directory
set NUMT = 16
set MEM = "32G"


# Path to the updated Python script
set PYTHON_SCRIPT = "/FULLPATH/TO/org_sep.py"

# Interactive sample selection
set s = ""

while ($s == "")
    echo ""
    echo "*********** Please enter a unique sample prefix for BAM separation, type ALL for all samples, or type exit to finish ***************"
    set s = "$<"
    if ($s == "exit") then
        continue
    else if ($s == "") then
        continue
    else if ($s == "ALL") then
        # Select all BAM files in the current directory
        set bam_files = `find ${BAM_DIR} -maxdepth 1 -name "*_outputAligned.sortedByCoord.out.bam"`
    else
        # Select BAM files matching a specific prefix
        set bam_files = `find ${BAM_DIR} -maxdepth 1 -name "${s}*_outputAligned.sortedByCoord.out.bam"`
        if (${#bam_files} == 0) then
            echo "No BAM files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through BAM files and submit jobs using sbatch
    foreach bam_file ($bam_files)
        set base = `basename ${bam_file} .bam`
        set output_base = `basename ${bam_file} _Caenorhabditis_elegans_miRNA-Seq_outputAligned.sortedByCoord.out.bam`
        set output_log = "${LOGS}/orgsep_${base}.txt"

        echo ""
        echo "********** Separating BAM by organism for sample $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              samtools sort -@ '${NUMT}' -n -O SAM '${bam_file}' > '${BAM_DIR}'/'${base}'.sam ; \
              python3 '${PYTHON_SCRIPT}' '${BAM_DIR}'/'${base}'.sam ; \
              samtools view -hbS '${BAM_DIR}'/'${output_base}'_ecoli.sam | samtools sort -@ '${NUMT}' -o '${BAM_DIR}'/'${output_base}'_ecoli.bam - ; \
              samtools index '${BAM_DIR}'/'${output_base}'_ecoli.bam ; \
              # Keep C. elegans SAM as SAM without conversion \
              samtools view -hbS '${BAM_DIR}'/'${output_base}'_ambiguous.sam | samtools sort -@ '${NUMT}' -o '${BAM_DIR}'/'${output_base}'_ambiguous.bam - ; \
              samtools index '${BAM_DIR}'/'${output_base}'_ambiguous.bam ; \
              rm -f '${BAM_DIR}'/'${output_base}'_ecoli.sam '${BAM_DIR}'/'${output_base}'_ambiguous.sam ;' \
              >! tmp_orgsep_${output_base}

        # Submit the job using sbatch and store the job ID for future dependency management
        sbatch -c ${NUMT} --mem=${MEM} --time=1-0 -o ${output_log} tmp_orgsep_${output_base} >&! ${BAM_DIR}/jobid_${output_base}

        # Wait briefly to avoid overwhelming the scheduler
        sleep 1
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_orgsep_*
echo "All BAM files processed and submitted successfully!"
