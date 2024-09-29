#!/bin/bash

#I use variables, for readability
#REF will be taken from Dan's directory
#vcf will be the file containing all our sequencing data (I prefer using the file without annotation, it's lighter this way; we will annotate later)
#output name your file
#bed contains the coordinates of all the genes (regions) of interest to be extracted
#sample_list containts the list of all the samples to be extracted (all cohorts together)


BASE_DIR=$1
REF=~/runs/go_lab/Reference/human_g1k_v37.fasta
vcf=$6
output_old=$BASE_DIR/$2.old.vcf
output=$BASE_DIR/$2.vcf
bed=$3
sample_list=$4
core=$5
mem=`echo "-Xmx"$((4*core))g`
java $mem -jar /lustre03/project/6004655/COMMUN/soft/lib/java/GATK/GenomeAnalysisTK-3.8/dist/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V $vcf -o $output_old -L $bed -sf $sample_list -env -nt $core --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES


#gatk command switches, explained:
# -T = which gatk algorithm to use
# -V = variant file (i.e. vcf); you can supply this one multiple times, e.g. -V file1.vcf -V file2.vcf etc
# -o = output. duh
# -L = intervals file; can be a bed file (chrom, start, stop), or a list of intervals: e.g. 1:1000  2:1238 2:1239 2:1240 etc (one per line)
# -sf = sample file; this is a list of samples (one per line) and MUST be written exactly as they appear in the vcf file header (on the line that starts with "#CHROM")
# -env - excllude non-variants; if you've removed samples with unique variants, this avoids having lines with no variant genotypes
