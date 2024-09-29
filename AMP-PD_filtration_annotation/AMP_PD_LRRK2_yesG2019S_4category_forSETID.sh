#!/bin/bash

module load nixpkgs/16.09
module load gcc/7.3 r/3.5.2

BASE_DIR=~/runs/sitkicem/JULY2024_yesG2019S_LRRK2


## FILTRATION 
plink --bfile ~/runs/senkkon/2022/TIM_TOM/AMP_PD/AMP_PD_FILTERED_ALL_CHR --extract range LRRK2.bed --keep covar_AMP_PD.txt --make-bed --out AMP_PD_all_variants

plink --bfile AMP_PD_all_variants --max-maf 0.01 --write-snplist --out AMP_PD_rare_variants

echo "rs34637584" >> AMP_PD_rare_variants.snplist

plink --bfile AMP_PD_all_variants --extract AMP_PD_rare_variants.snplist --make-bed --out AMP_PD

plink --bfile ~/runs/senkkon/2022/TIM_TOM/AMP_PD/AMP_PD_FILTERED_ALL_CHR --extract range LRRK2.bed --keep covar_AMP_PD.txt --max-maf 0.01 --make-bed --out AMP_PD


## ANNOTATION

awk 'FNR==NR { a[$2]; next } ($6 in a) { print $0 }' AMP_PD.bim ~/runs/sitkicem/annovar/humandb/hg38_avsnp150.txt > annotation/AMP_PD.avinput

perl ~/runs/sitkicem/annovar/table_annovar.pl annotation/AMP_PD.avinput ~/runs/sitkicem/annovar/humandb/ --buildver hg38 --out annotation/AMP_PD.avinput_recode_convert.annovar --remove --protocol refGene,ljb26_all,dbnsfp41c,avsnp150 --operation g,f,f,f --nastring .

## CATEGORIZATION

# awk 'BEGIN {FS=OFS="\t"} {if (NR==1 || ($7 == "LRRK2" && ($9 ~ /stop/ || $9 ~ /nonsyn/ || $9 ~ /frame/ || ($9 ~ /intronic_splicing/ && $9 ~ /[+-][1-2]/)))) print $79}' AMP_PD.avinput_recode_convert.annovar.hg38_multianno.txt > AMP_PD_all_funct.txt

# NON-SYNONYMOUS VARIANTS
awk -F "\t" '{ if (NR==1 ||$9 ~ /nonsyn/ ) print $79}' annotation/AMP_PD.avinput_recode_convert.annovar.hg38_multianno.txt > annotation/AMP_PD_nonsyn.txt

# LOF VARIANTS
awk -F "\t" '{ if (NR==1 || $9 ~ /stop/ || $9 ~ /frame/ || ($9 ~ /intronic_splicing/ && $9 ~ /[+-][1-2]/) ) print $79}' annotation/AMP_PD.avinput_recode_convert.annovar.hg38_multianno.txt > annotation/AMP_PD_LOF.txt

# ENCODE VARIANTS
awk -F "\t" '{ if (NR==1 || $6 == "exonic") print $79}' annotation/AMP_PD.avinput_recode_convert.annovar.hg38_multianno.txt > annotation/AMP_PD_ENCODE.txt

# CADD VARIANTS
awk -F "\t" '{ if (NR==1 || $31 >= 20) print $79}' annotation/AMP_PD.avinput_recode_convert.annovar.hg38_multianno.txt > annotation/AMP_PD_CADD.txt
