#!/bin/tcsh
###############################################################################
# STAR Alignment Pipeline
# Author: Noa Gilad
#
# Description:
#   Interactive script for aligning paired-end FASTQ files with STAR.
#   Submits each sample as an sbatch job with dependency on a previous job.
#
# Usage:
#   The script will prompt for:
#       - STAR genome index directory
#       - Job directory containing previous job IDs
#       - Job ID prefix (used for dependency)
#       - FASTQ file suffix (e.g., "1" or "_R1.fastq")
#       - Full suffix for removing basename (e.g., "_R1.fastq")
#       - Full suffix for R2 reads
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - STAR module available
#
# User-configurable paths:
#   LOGS        → Directory for SLURM log files
###############################################################################

set USAGE = "Usage: STAR: an interactive script, will ask for the needed parameters"
set tooMany = $1

if ($tooMany != "") then
    echo $USAGE
    exit 1;
endif

### PARAMETERS ###
# Set working directory as the running folder
set DIR = `pwd`
set OUT_DIR = "./STAR_out"
set LOGS = "/FULLPATH/TO/LOGS"  # Specified SLURM log directory
set NUMT = 12
set MEM = "128G"

echo "*********** Please enter index directory ***************"
set INDEX = "$<"

echo "*********** Please enter job dir directory ***************"
set JOB_DIR = "$<"

echo "*********** Please enter job id prefix ***************"
set JOB_ID = "$<"

echo "*********** Please enter file suffix (for paired choose 1) ***************"
set FILE_SUFFIX = "$<"

echo "*********** Please enter file full suffix for basename removal ***************"
set FULL_SUFFIX = "$<"

echo "*********** Please enter file full suffix of R2 ***************"
set FULL_SUFFIX_R2 = "$<"

mkdir ${OUT_DIR}

# Interactive sample selection
set s = ""

while ($s == "")
    echo ""
    echo "*********** Please enter a unique sample prefix for fastq, type ALL for all samples, or type exit to finish ***************"
    set s = "$<"
    if ($s == "exit") then
        continue
    else if ($s == "") then
        continue
    else if ($s == "ALL") then
        # Select all fastq files in the current directory
        set fastq_files = `find ${DIR} -maxdepth 1 -name "*${FILE_SUFFIX}"`
    else
        # Select fastq files matching a specific prefix
        set fastq_files = `find ${DIR} -maxdepth 1 -name "${s}*${FILE_SUFFIX}"`
        if (${#fastq_files} == 0) then
            echo "No fastq files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through fastq files and submit jobs using sbatch
    foreach fastq_file ($fastq_files)
        set base = `basename ${fastq_file} ${FULL_SUFFIX}`
        set output_log = "${LOGS}/star_${base}.txt"
        set output_file = "${OUT_DIR}/${base}"

        echo ""
        echo "********** aligning $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              module load star ; \
              STAR --runMode alignReads --genomeDir '${INDEX}' --outSAMtype BAM SortedByCoordinate --readFilesIn '${DIR}'/'${base}''${FULL_SUFFIX}' '${DIR}'/'${base}''${FULL_SUFFIX_R2}'  --runThreadN 12 --outFileNamePrefix '${output_file}' ;' \
              >! tmp_star_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${JOB_DIR}/${JOB_ID}${base} | cut -d ' ' -f 4`

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=2-0 -o ${output_log} -d afterok:$alid tmp_star_${base}) >&! ${DIR}/jobid_star_${base}


        # Wait briefly to avoid overwhelming the scheduler
        sleep 2
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_star_*
echo "All fastq files processed and submitted successfully!"
