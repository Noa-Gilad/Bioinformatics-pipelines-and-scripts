#!/bin/tcsh
###############################################################################
# Smart-seq3 Trimming Pipeline (Cutadapt-based)
# Author: Noa Gilad
# Description:
#   Interactive trimming script for Smart-seq3 paired-end FASTQ files.
#   Supports polyA trimming, adapter trimming, and ligation-based trimming.
#   Submits each sample as a SLURM job using sbatch.
#
# Usage:
#   Run the script inside a directory containing FASTQ files named:
#       SAMPLE_R1_001.fastq
#       SAMPLE_R2_001.fastq
#   Or change the suffix in this script.
#
#   The script will prompt the user for:
#       - A specific sample prefix
#       - OR "ALL" to process all samples
#
# Requirements:
#   - SLURM scheduler available (`sbatch`)
#   - cutadapt installed and loadable via module OR conda
#
# Paths requiring user configuration:
#   LOGS       → Path to desired log directory
#   ADAPTERS   → Path to FASTA file containing adapter sequences
###############################################################################

set USAGE = "Usage: trimming: an interactive script, will ask for the needed parameters"
set tooMany = $1

if ($tooMany != "") then
    echo $USAGE
    exit 1;
endif

### PARAMETERS ###
# Set working directory as the running folder
set DIR = `pwd`
set OUT_DIR = "./trimmed"
set LOGS = "/FULLPATH/TO/LOGS"  # Specified SLURM log directory
set NUMT = 16
set MEM = "128G"
set ADAPTERS = "/FULLPATH/TO/polyN/ADAPTERS.fa"
set POLY_DIR = "./trimmed/poly_trimmed/"
set LIG_DIR = "./trimmed/poly_trimmed/ligation"
set NO_LIG_DIR = "./trimmed/poly_trimmed/no_ligation"

mkdir -p ${LIG_DIR}
mkdir -p ${NO_LIG_DIR}

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
        set fastq_files = `find ${DIR} -maxdepth 1 -name "*_R1_001.fastq"`
    else
        # Select fastq files matching a specific prefix
        set fastq_files = `find ${DIR} -maxdepth 1 -name "${s}*_R1_001.fastq"`
        if (${#fastq_files} == 0) then
            echo "No fastq files found for the prefix '$s'. Please try again."
            set s = ""
            continue
        endif
    endif

    # Loop through fastq files and submit jobs using sbatch
    foreach fastq_file ($fastq_files)
        set base = `basename ${fastq_file} _R1_001.fastq`
        set output_log = "${LOGS}/cutadapt_${base}.txt"

        echo ""
        echo "********** trimming $base ***********"

        # Create a temporary script for sbatch submission
        echo '#\!/bin/tcsh \
              module load py-cutadapt ; \
              cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -m 15 -o '${OUT_DIR}'/'${base}'_trimmed_1.fastq -p '${OUT_DIR}'/'${base}'_trimmed_2.fastq '${DIR}'/'${base}'_R1_001.fastq '${DIR}'/'${base}'_R2_001.fastq ; \
              cutadapt -a file:'${ADAPTERS}' -A file:'${ADAPTERS}' -m 15 -o '${POLY_DIR}'/'${base}'_trimmed_poly_1.fastq -p '${POLY_DIR}'/'${base}'_trimmed_poly_2.fastq '${OUT_DIR}'/'${base}'_trimmed_1.fastq '${OUT_DIR}'/'${base}'_trimmed_2.fastq ; \
              cutadapt -j 16 -g "^ATTGCGCAATG" -O 8 -m 15 --untrimmed-output '${NO_LIG_DIR}'/'${base}'_nolig_R1.fastq --untrimmed-paired-output '${NO_LIG_DIR}'/'${base}'_nolig_R2.fastq -o '${LIG_DIR}'/'${base}'_trimmed_lig_R1.fastq -p '${LIG_DIR}'/'${base}'_trimmed_lig_R2.fastq '${POLY_DIR}'/'${base}'_trimmed_poly_1.fastq '${POLY_DIR}'/'${base}'_trimmed_poly_2.fastq ;' \
              >! tmp_trimming_${base}

        # Submit the job using sbatch with a dependency on the previous job ID
        (sbatch -c ${NUMT} --mem=${MEM} --time=1-0 -o ${output_log} tmp_trimming_${base}) >&! ${DIR}/jobid_trimming_${base}


        # Wait briefly to avoid overwhelming the scheduler
        sleep 2
    end

    # Prepare for next sample or exit the loop
    if ($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_trimming_*
echo "All fastq files processed and submitted successfully!"
