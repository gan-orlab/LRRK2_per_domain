#!/bin/bash

module load StdEnv/2020
module load gatk/4.2.5.0

BASE_DIR=$1
vcf=$2
GQ=$3
DP=$4
AN=$5
ANratio=$6
bed_file=$7
output=LRRK2


hg38ref=/lustre03/project/6004655/COMMUN/runs/senkkon/UKBB_DP/GRCh38_full_analysis_set_plus_decoy_hla.fa
# SET INTERVAL FOR GENE (LRRK2 38 MINUTES)
gatk SelectVariants \
  -R $hg38ref \
  -V ${vcf}.vcf.gz \
  -O ${output}.vcf.gz \
  -L ${bed_file} \
  --java-options "-Xmx7g" &&

#REMOVES GENOTYPE CALLS PER SAMPLE THAT ARE LOWER THAN GQ20 (LRRK2 43 MINUTES)
# gatk VariantFiltration \
  -R $hg38ref \
  -V ${output}.vcf.gz \
  -O ${output}_GQ${GQ}.vcf.gz \
  --java-options "-Xmx7g" \
  --genotype-filter-expression "GQ < ${GQ}" \
  --genotype-filter-name "GQ${GQ}" \
  --set-filtered-genotype-to-no-call \
  --missing-values-evaluate-as-failing &&

#REMOVES GENOTYPE CALLS PER SAMPLE THAT ARE LOWER THAN DP10 (LRRK2 45 MINUTES)
# gatk VariantFiltration \
  -R $hg38ref \
  -V ${output}_GQ${GQ}.vcf.gz \
  -O ${output}_GQ${GQ}_DP${DP}.vcf.gz \
  --java-options "-Xmx12g" \
  --genotype-filter-expression "DP < ${DP}" \
  --genotype-filter-name "DP${DP}" \
  --set-filtered-genotype-to-no-call \
  --missing-values-evaluate-as-failing &&

# FLAGS POSITIONS AS "TO REMOVE" IF LESS THAN THE ALLOWABLE NUMBER OF ALLELES ARE PRESENT AFTER ABOVE FILTERING
gatk VariantFiltration \
  -R $hg38ref \
  -V ${output}_GQ${GQ}_DP${DP}.vcf.gz \
  -O ${output}_GQ${GQ}_DP${DP}_MISS${ANratio}.vcf.gz \
  --java-options "-Xmx12g" \
  --filter-expression "AN < $AN" \
  --filter-name "MISS" \
  --missing-values-evaluate-as-failing &&

# REMOVE THE FLAGGED POSITIONS FROM THE VCF AND OUTPUT A CLEANED VERSION
gatk SelectVariants \
  -R $hg38ref \
  -V ${output}_GQ${GQ}_DP${DP}_MISS${ANratio}.vcf.gz \
  -O ${output}_GQ${GQ}_DP${DP}_MISS${ANratio}_filtered.vcf.gz \
  --java-options "-Xmx12g" \
  --exclude-filtered \
  --exclude-non-variants &&

#CLEAN UP INTERMEDIATE FILES
rm ${output}.vcf.gz ${output}.vcf.gz.tbi
rm ${output}_GQ${GQ}.vcf.gz ${output}_GQ${GQ}.vcf.gz.tbi
rm ${output}_GQ${GQ}_DP${DP}.vcf.gz ${output}_GQ${GQ}_DP${DP}.vcf.gz.tbi
rm ${output}_GQ${GQ}_DP${DP}_MISS${ANratio}.vcf.gz ${output}_GQ${GQ}_DP${DP}_MISS${ANratio}.vcf.gz.tbi &&

#ALL DONE
echo "filtration of variants in ${vcf} is a great success !!"
