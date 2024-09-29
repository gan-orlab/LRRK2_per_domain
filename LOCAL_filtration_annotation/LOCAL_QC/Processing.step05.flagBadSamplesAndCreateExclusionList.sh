#!/bin/bash


###Same as script for step05 but I decided to go with a threshold of 50% for missing genotyping and added a line to produce the exclusion list for this threshold
 
vcf=$1
out=${2:-$(basename $vcf).sample2fractionNoCALLs}
out2=$(basename $vcf).50PercentShitSamplesToExclude
out3=$(basename $vcf).10PercentShitSamplesToExclude
if [[ ! -s $vcf ]]; then echo "WTF COWBOY YOUR INPUT FILE (NAME=$vcf) IS EMPTY OR DOES NOT EXIST"; exit 42; fi

awk -F"\t" 'BEGIN{FS=OFS="\t"}{if ($0~/^#/) {if ($1=="#CHROM") {for (i=10; i<=NF; i++) samples[i]=$i}; next}; total++; for (i=10; i<=NF; i++) if ($i~/^\.\/\./) c[samples[i]]++; }END{for (i in c) print i,c[i]/total}' $vcf |sort -k2n > $out


cat $out | awk '{ if ($2>0.5) {print $1}}'  > $out2
cat $out | awk -v cut=$3 '{ if ($2>0.1) {print $1}}'  > $out3
