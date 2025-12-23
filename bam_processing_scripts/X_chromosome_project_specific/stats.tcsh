#!/bin/tcsh
###############################################################################
# stats.tcsh
# Author: Noa Gilad
#
# Description:
#   Interactive wrapper to compute per-sample statistics for the
#   X-chromosome project, submitting one SLURM job per input file.
#
# Usage:
#   Run inside a directory containing the relevant input files and follow
#   the prompts for a sample prefix or "ALL".
#
# Requirements:
#   - SLURM scheduler (`sbatch`)
#   - A gene-biotype annotation table (BIO_FILE) for downstream parsing.
#
# Paths requiring user configuration:
#   LOGS     → directory where SLURM log files will be written
#   BIO_FILE → full path to gene_biotype.sorted.txt (or equivalent)
###############################################################################

set USAGE = "Usage: stats: an interactive script, will ask for the needed parameters"
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
set BIO_FILE = "/FULLPATH/TO/gene_biotype.sorted.txt"

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

    # Loop through BAM files and run htseq-count
    foreach bam_file ($bam_files)
        set base = `basename ${bam_file} _Caenorhabditis_elegans_miRNA-Seq_outputAligned.sortedByCoord.out.bam`
        set output_file = "${SAM_DIR}/${base}_stats.txt"
        set output_log = "${LOGS}/run_stats_${base}.txt"

        echo ""
        echo "********** Analysing stats on sample $base ***********"

        echo '#\!/bin/tcsh \
              echo "Total and Mapped Reads Per BAM File" >! '${output_file}' \
              set total = `samtools view -c '${bam_file}'` \
              set mapped = `samtools view -c -F 4 '${bam_file}'` \
              set mapped_ecoli = `samtools view -c '${SAM_DIR}'/'${base}'_ecoli.bam` \
              set ambi_org = `samtools view -c '${SAM_DIR}'/'${base}'_ambiguous.bam` \
              set cel_org = `samtools view -c '${SAM_DIR}'/'${base}'_celegans.sam` \
              set final = `samtools view -c '${SAM_DIR}'/'${base}'_cel.bam` \
\
              set ecoli_per = `echo "scale=2; $mapped_ecoli*100/$mapped" | bc` \
              set ambi_org_per = `echo "scale=2; $ambi_org*100/$mapped" | bc` \
              set cel_org_per = `echo "scale=2; $cel_org*100/$mapped" | bc` \
              set final_from_mapped = `echo "scale=2; $final*100/$mapped" | bc` \
              set final_from_total = `echo "scale=2; $final*100/$total" | bc` \
\
              echo "\\nOriginal bam total reads\\t$total" >> '${output_file}'  \
              echo "\\nOriginal bam mapped reads\\t$mapped" >> '${output_file}'  \
\
              echo  "\\nMapped reads to e.coli\\t$mapped_ecoli" >> '${output_file}' \
              echo  "\\nPercentage of mapped reads to e.coli from all mapped reads\\t$ecoli_per%" >> '${output_file}' \
\
              echo  "\\nMapped reads to ambiguous organisms\\t$ambi_org" >> '${output_file}' \
              echo  "\\nPercentage of reads mapped both to e.coli and to c.elegans from all mapped reads\\t$ambi_org_per%" >> '${output_file}' \
\
              echo  "\\nMapped reads to celegans\\t$cel_org" >> '${output_file}' \
              echo  "\\nPercentage of reads mapped only to c.elegans from all mapped reads\\t$cel_org_per%" >> '${output_file}' \
\
              echo "\\nMapped reads to c.elegans only and has the highest score\\t$final" >> '${output_file}' \
              echo  "\\nPercentage of reads mapped to c.elegans only and has the highest score from all mapped reads\\t$final_from_mapped%" >> '${output_file}' \
              echo  "\\nPercentage of reads mapped to c.elegans only and has the highest score from total reads\\t$final_from_total%" >> '${output_file}' \
\
              echo "\\nTotal Reads Mapped Per Biotype" >> '${output_file}' \
\
              echo "\\nBiotype\\tTotal_Reads" >> '${output_file}' \
              sort -k1,1 '${SAM_DIR}'/'${base}'_counts.txt >! '${SAM_DIR}'/'${base}'_counts.sorted.txt \
\
              join -1 1 -2 1 '${BIO_FILE}' '${SAM_DIR}'/'${base}'_counts.sorted.txt | awk '\''{sum[$2] += $3} END {for(b in sum) print b", "sum[b]}'\'' >> '${output_file}' \
\
              echo "\\nGenes With More Than 2 Reads Per Biotype" >> '${output_file}' \
              echo "\\nBiotype\tGenes_With_More_Than_2_Reads" >> '${output_file}' \
\
              join -1 1 -2 1 '${BIO_FILE}' '${SAM_DIR}'/'${base}'_counts.sorted.txt | awk '\''{if($3 > 2) count[$2]++} END {for(b in count) print b", "count[b]}'\'' >> '${output_file}' \
\
              ' >! tmp_stats_${base}

        set alid = `less ${SAM_DIR}/jobid_htseq_${base} | cut -d ' ' -f 4`
        (sbatch -c ${NUMT} --mem=${MEM} --time=2-0 --exclude=glacier-20 -o ${output_log} -d afterok:$alid tmp_stats_${base}) >&! ${SAM_DIR}/jobid_stats_${base}
        sleep 1
    end

    # Prepare for next sample
    if($s != "ALL") then
        set s = ""
    endif
end

wait
rm -f tmp_stats_*
echo "stats completed!"
