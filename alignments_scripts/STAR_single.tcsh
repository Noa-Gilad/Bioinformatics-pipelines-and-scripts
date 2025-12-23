#!/bin/tcsh
###############################################################################
# STAR_single.tcsh
# Author: Noa Gilad
#
# Description:
#   Interactive STAR alignment wrapper for single-end FASTQ files.
#   For each FASTQ in the current directory, submits a SLURM job that runs:
#     STAR --runMode alignReads --genomeDir INDEX --outSAMtype BAM SortedByCoordinate
#
# Usage:
#   Run inside a directory containing *.fastq files.
#   The script will prompt for a sample prefix or "ALL" to process all files.
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - STAR installed and available (e.g. via `module load star`)
#
# Paths requiring user configuration:
#   LOGS  → directory where SLURM log files will be written
#   INDEX → STAR genomeDir path
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
set LOGS = "/FULLPATH/TO/LOGS"        # Specified SLURM log directory
set INDEX = "/FULLPATH/TO/STAR_INDEX" # STAR genomeDir path
set JOB_DIR = "/FULLPATH/TO/JOB_DIR" # The directory in which the former jobs are in
set NUMT = 12
set MEM = "128G"

mkdir -p ${OUT_DIR}

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
        set fastq_files = `find ${DIR} -maxdepth 1 -name "*.fastq"`
    else
        # Select fastq files matching a specific prefix
        set fastq_files = `find ${DIR} -maxdepth 1 -name "${s}*.fastq"`
        if (${#fastq_files} == 0) then
            echo "No fastq files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through fastq files and submit jobs using sbatch
    foreach fastq_file ($fastq_files)
        set base = `basename ${fastq_file} _trimmed_poly.fastq`
        set output_log = "${LOGS}/star_${base}.txt"
        set output_file = "${OUT_DIR}/${base}"

        echo ""
        echo "********** aligning $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              module load star ; \
              STAR --runMode alignReads --genomeDir '${INDEX}' --outSAMtype BAM SortedByCoordinate --readFilesIn '${DIR}'/'${base}'_trimmed_poly.fastq --runThreadN 12 --outFileNamePrefix '${output_file}' ;' \
              >! tmp_star_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${JOB_DIR}/jobid_trimming_${base} | cut -d ' ' -f 4`

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=2-0  --exclude=glacier-33 -o ${output_log} -d afterok:$alid tmp_star_${base}) >&! ${DIR}/jobid_star_${base}


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
