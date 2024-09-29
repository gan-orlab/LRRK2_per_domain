#!/bin/bash

#I use variables, for readability
#####LOOPING NOT DONE YET, WORK ON IT WITH DAN
###general code; can be put in a loop if you have a file with list of genes + specify the cohorts: for gene in selectedgenes.list; do for cohort in FC NY ISR; do for DP in 15 30 50; do bash script.sh $gene $cohort $DP; done; done; done 

###Khair's version: can be put in a loop if you have a file with list of genes: for gene in PARK2 PARK7 PINK1 VPS13C; do for DP in 15 30 50; do bash script.sh $gene RBD $DP; done; done




read vcf cohort_name DP DIR core<<< $@
if [[ ! -s $vcf ]]; then echo "ERROR: input vcf (1st arg) not specified, empty or does not exist"; exit 42; fi
if [[ -z $DP || ! $DP -gt 0 ]]; then echo "ERROR: cohort name (3rd arg) not specified"; exit 42; fi
if [[ -z $cohort_name ]]; then echo "ERROR: cohort name (3rd arg) not specified"; exit 42; fi

REF=~/runs/go_lab/Reference/human_g1k_v37.fasta

sample_list=$DIR/$cohort_name.samples.list
output=$cohort_name.DP$DP.all.genes.vcf
mem=`echo "-Xmx"$((core*4))g`
java $mem -jar /lustre03/project/6004655/COMMUN/soft/lib/java/GATK/GenomeAnalysisTK-3.8/dist/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V $vcf -o $output -sf $sample_list -env --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES -nt $core

#for cohort in FC NY ISR; do bash ../Scripts/Processing.step07.selectByCohorts.sh PD_familial_genes_except_GBA_AllSamples_GF25_annotated_GQ30_DP30_MISS10_filtered_cleaned.vcf $cohort 30; done
