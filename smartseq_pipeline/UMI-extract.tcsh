#!/bin/tcsh
###############################################################################
# UMI Extraction Pipeline (umi-tools)
# Author: Noa Gilad
#
# Description:
#   Interactive script for extracting UMIs from Smart-seq3 ligation-trimmed reads.
#   Works on paired FASTQ files ending in:
#       SAMPLE_trimmed_lig_R1.fastq
#       SAMPLE_trimmed_lig_R2.fastq
#
# Usage:
#   Run the script in a directory containing the FASTQ files above.
#   The script will prompt for:
#       - A specific sample prefix
#       - OR "ALL" to process all samples
#
# Requirements:
#   - SLURM scheduler available (`sbatch`)
#   - umi_tools installed via conda
#
# Paths requiring user configuration:
#   LOGS     → Directory where SLURM logs should be written
#   JOB_DIR  → Directory where jobid_trimming_* files are stored
###############################################################################

set USAGE = "Usage: umi extraction: an interactive script, will ask for the needed parameters"
set tooMany = $1

if ($tooMany != "") then
    echo $USAGE
    exit 1;
endif

### PARAMETERS ###
# Set working directory as the running folder
set DIR = `pwd`
set OUT_DIR = "./umi_extracted"
set LOGS = "/FULLPATH/TO/LOGS"  # Specified SLURM log directory
set NUMT = 16
set MEM = "32G"
set JOB_DIR = "/FULLPATH/TO/jobid_trimming_DIRECTORY"

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
        set fastq_files = `find ${DIR} -maxdepth 1 -name "*_R1.fastq"`
    else
        # Select fastq files matching a specific prefix
        set fastq_files = `find ${DIR} -maxdepth 1 -name "${s}*_R1.fastq"`
        if (${#fastq_files} == 0) then
            echo "No fastq files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through fastq files and submit jobs using sbatch
    foreach fastq_file ($fastq_files)
        set base = `basename ${fastq_file} _trimmed_lig_R1.fastq`
        set output_log = "${LOGS}/umi_extraction_${base}.txt"

        echo ""
        echo "********** extracting UMI from $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              module load spack miniconda3 ; \
              source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.csh ; \
              conda activate umi-tools-env ; \
              umi_tools extract -I '${DIR}'/'${base}'_trimmed_lig_R1.fastq --bc-pattern=NNNNNNNN --read2-in='${DIR}'/'${base}'_trimmed_lig_R2.fastq --stdout='${OUT_DIR}'/'${base}'_UMI_extracted_R1.fastq --read2-out='${OUT_DIR}'/'${base}'_UMI_extracted_R2.fastq ;' \
              >! tmp_umi_extract_${base}

        # Extract the previous job ID for dependency (without using $d)
        set alid = `less ${JOB_DIR}/jobid_trimming_${base} | cut -d ' ' -f 4`

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=1-0 -o ${output_log} -d afterok:$alid tmp_umi_extract_${base}) >&! ${DIR}/jobid_umi_extraction_${base}


        # Wait briefly to avoid overwhelming the scheduler
        sleep 2
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_umi_extract_*
echo "All fastq files processed and submitted successfully!"
