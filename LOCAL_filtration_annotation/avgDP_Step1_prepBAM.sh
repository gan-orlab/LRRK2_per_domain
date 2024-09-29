#!/bin/bash

# # Define the base directory where the BAM files are stored
# base_dir="/lustre03/project/6004655/COMMUN/pipeline_results/mip/PD/PD_MIP71/PD_MIP71"

# # Define the sequencing information and processing steps subdirectories
# sequencing_info="Illumina_HiSeq_Paired-IC-MIP-2023_06_28/GATK_BWA.v37"

all_bam_file="/lustre03/project/6004655/COMMUN/runs/go_lab/GRCh37/mips/unfiltered/PD/all_bam_file_up-to-MIP103.txt"

# Define the list of cohorts and their corresponding sample files
declare -A cohorts
cohorts=( ["FRENCH"]="cohort_FRENCH/FRENCH.samples.list" ["RUSSIA"]="cohort_RUSSIA/RUSSIA.samples.list" ["ISRAEL"]="cohort_ISRAEL/ISRAEL.samples.list" ["USA"]="cohort_USA/USA.samples.list" )

# Iterate over each cohort and process the sample file
for cohort in "${!cohorts[@]}"; do
    samples_file=${cohorts[$cohort]}
    
    echo "Processing cohort: $cohort"
    
    # Create or empty the cohort BAM list file
    > "${cohort}_BAM.list"
    
    valid_count=0

    # Read the sample identifiers from the file and search for BAM file paths
    while IFS= read -r sample_id; do
        bam_file=$(grep "$sample_id" "$all_bam_file")
        
        if [ -n "$bam_file" ]; then
            echo "$bam_file" >> "${cohort}_BAM.list"
            valid_count=$((valid_count + 1))
        fi
    done < "$samples_file"
    
    echo "Found $valid_count valid BAM files for cohort: $cohort"
    echo ""
done
