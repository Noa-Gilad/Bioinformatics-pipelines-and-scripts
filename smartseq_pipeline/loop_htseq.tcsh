#!/bin/tcsh
###############################################################################
# loop_htseq.tcsh
# Author: Noa Gilad
#
# Description:
#   Interactive wrapper for running htseq-count on multiple BAM files.
#   Supports dependency handling based on previous pipeline steps (STAR/UMI).
#
# Usage:
#   Run inside directory containing BAM files.
#   The script will ask interactively for:
#       - job directory
#       - GTF file
#       - job ID prefix
#       - BAM suffixes
#       - strandedness mode
#       - sample prefix or "ALL"
#
# Requirements:
#   - SLURM scheduler (sbatch)
#   - HTSeq module available (py-htseq)
#
###############################################################################

set USAGE = "Usage: loop_htseq: an interactive script, will ask for the needed parameters"
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
set MEM = "64g"
set HTSEQ_MODULE = "py-htseq"

echo "*********** Please enter job DIR ***************"
set JOB_DIR = "$<"

echo "*********** Please enter GTF path ***************"
set GTF_FILE = "$<"

echo "*********** Please enter job id prefix ***************"
set JOB_ID = "$<"

echo "*********** Please enter file suffix ***************"
set FILE_SUFFIX = "$<"

echo "*********** Please enter file full suffix for basename removal ***************"
set FULL_SUFFIX = "$<"

echo "**** Please enter stranded parameter - yes for stranded, no for unstranded, reverse for stranded with reverse ****"
set STRANDED = "$<"


# Interactive sample selection
set s = ""

while($s == "")
    echo ""
    echo "*********** Please enter a unique sample prefix for BAM, type ALL for all samples, or type exit to finish ***************"
    set s = "$<"
    if ($s == "exit") then
        continue
    else if ($s == "") then
        continue
    else if ($s == "ALL") then
        # Select all BAM files in the current directory
        set bam_files = `find ${SAM_DIR} -maxdepth 1 -name "*${FILE_SUFFIX}"`
    else
        # Select BAM files matching a specific prefix
        set bam_files = `find ${SAM_DIR} -maxdepth 1 -name "${s}*${FILE_SUFFIX}"`
        if (${#bam_files} == 0) then
            echo "No BAM files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through BAM files and run htseq-count
    foreach bam_file ($bam_files)
        set base = `basename ${bam_file} ${FULL_SUFFIX}`
        set output_file = "${SAM_DIR}/${base}_counts.txt"
        set output_log = "${LOGS}/run_htseq_${base}.txt"

        echo ""
        echo "********** Running HTSeq on sample $base ***********"

        echo "#\!/bin/tcsh \
              module load ${HTSEQ_MODULE} \
              htseq-count -f bam -r pos -s '${STRANDED}' '${SAM_DIR}'/'${base}''${FULL_SUFFIX}' '${GTF_FILE}' >! '${output_file}' \
              " >! tmp_htseq_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${JOB_DIR}/${JOB_ID}${base} | cut -d ' ' -f 4`

        (sbatch -c ${NUMT} --mem=${MEM} --time=2-0 --exclude=glacier-04 -o ${output_log} -d afterok:$alid tmp_htseq_${base}) >&! ${SAM_DIR}/jobid_htseq_${base}
        sleep 2
    end

    # Prepare for next sample
    if($s != "ALL") then
        set s = ""
    endif
end

wait
echo "HTSeq-count completed!"
