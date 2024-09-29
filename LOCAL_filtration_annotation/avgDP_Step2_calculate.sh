#!/bin/bash

lrrk2_bed="LRRK2_hg19.bed"


# Process each cohort
calculate_total_depth() {
    local bam_file=$1
    local lrrk2_bed=$2
    
    # Use samtools depth to calculate depth
    samtools depth -b "$lrrk2_bed" "$bam_file" | awk '{sum+=$3; count++} END {print sum, count}'
}


# Process each cohort
for cohort in "USA" "ISRAEL" "FRENCH"; do
    bam_list="${cohort}_BAM.list"
    
    echo "Processing cohort: $cohort"
    
    total_depth=0
    total_count=0
    sample_count=$(wc -l < "$bam_list")

    total_files=$(wc -l < "$bam_list")
    processed_files=0

    # Calculate the average depth for each BAM file
    while IFS= read -r bam_file; do
        processed_files=$((processed_files + 1))
        
        # Calculate progress percentage every 50 iterations
        if [ $((processed_files % 50)) -eq 0 ]; then
            echo -ne "Calculating depth for $bam_file ($processed_files/$total_files) \r"
        fi

        read bam_total_depth bam_count < <(calculate_total_depth "$bam_file" "$lrrk2_bed")
        total_depth=$(awk "BEGIN {print $total_depth + $bam_total_depth}")
        total_count=$((total_count + bam_count))
        
        # Print the depth of the current BAM file
        echo "Depth for $bam_file: $(awk "BEGIN {print $bam_total_depth / $bam_count}")"
    done < "$bam_list"
    
    # Calculate the average depth across all samples in the cohort
    if [ $total_count -gt 0 ]; then
        average_depth=$(awk "BEGIN {print $total_depth / $total_count}")
    else
        average_depth=0
    fi
    
    echo "Cohort $cohort average depth: $average_depth"
    echo ""
    
done | tee average_depth_test3_slow.txt

# TESTING

samtools depth -b "LRRK2_hg19.bed" "/lustre03/project/6004655/COMMUN/pipeline_results/mip/PD/PD_MIP97/PD_MIP97/S40314/Illumina_HiSeq_Paired-IC-MIP-2022_05_25/GATK_BWA.v37/MIP.S40314.clean.bam" | awk '{sum+=$3} END {if (NR > 0) print "Average =", sum/NR; else print "Average = 0"}'
