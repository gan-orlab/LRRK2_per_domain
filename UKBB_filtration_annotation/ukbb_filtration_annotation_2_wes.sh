#!/bin/bash

BASE_DIR="~/runs/sitkicem/OCT_UKBB"

# STEP 0: Indexing the VCF file via tabix
srun -c 1 --mem=10g -t 0:15:0 bash tabix.sh WES_UKB_LRRK2.vcf.gz
# STEP 1: Limiting the VCF file by LRRK2 boundaries 40618781	40763172 and filtering for GQ20, DP10, 5% missingness
srun -c 1 --mem=20g -t 3:0:0 bash GATK_interval_and_filtration.sh ~/runs/sitkicem/OCT_UKBB WES_UKB_LRRK2 20 10 883500 95 LRRK2.bed 2>&1 | tee 26_OCT_GATK_MISS_log.txt &&
# STEP 1.5: Filter multi-allelic variants
bcftools view --max-alleles 2 -o LRRK2_GQ20_DP10_MISS95_filtered_noMultiAllelic.vcf.gz LRRK2_GQ20_DP10_MISS95_filtered.vcf.gz
# STEP 2: Convert VCF output after GATK into PLINK format
srun -c 4 --mem=15g -t 0:30:0 plink --vcf LRRK2_GQ20_DP10_MISS95_filtered_noMultiAllelic.vcf.gz --vcf-half-call m --make-bed --out LRRK2_after_GATK --output-chr M --max-maf 0.01 &&
# STEP 3: Convert VCF to ANNOVAR format
srun -c 1 --mem=20g -t 3:0:0 perl ~/runs/sitkicem/annovar/convert2annovar.pl --format vcf4 LRRK2_GQ20_DP10_MISS95_filtered_noMultiAllelic.vcf.gz --allsample --withfreq --outfile LRRK2_recode_convert &&
# STEP 4: Annotate SNPs with MAF3 using ANNOVAR
srun -c 1 --mem=20g -t 3:0:0 perl ~/runs/sitkicem/annovar/table_annovar.pl LRRK2_recode_convert ~/runs/sitkicem/annovar/humandb/ --buildver hg38 --out LRRK2_recode_convert.annovar --remove --protocol refGene,ljb26_all,dbnsfp41c --operation g,f,f --nastring .

# mkdir -p annotated
# mkdir -p annotated/all_coding
# mkdir -p annotated/CADD
# mkdir -p annotated/all_rare
# mkdir -p annotated/all_functional

# mkdir -p annotated_coding

# # # CODING NON-SYNONYMOUS VARIANTS
awk 'BEGIN {FS=OFS="\t"} $6 ~ /exonic/ && $9 ~ /nonsyn/ {print $1,$2,$3,$4,$5,$7}' LRRK2_recode_convert.annovar.hg38_multianno.txt > annotated_coding/LRRK2_all_coding.txt

# # Prepare SETID
awk 'BEGIN {FS=OFS="\t"} NR > 1 && $3 != "." {
                if ($3 >= 40225153 && $3 <= 40277914) domain="ARM";
                else if ($3 >= 40277915 && $3 <= 40284043) domain="ANK";
                else if ($3 >= 40294854 && $3 <= 40305881) domain="LRR";
                else if ($3 >= 40308486 && $3 <= 40322084) domain="ROC-COR";
                else if ($3 >= 40323249 && $3 <= 40351571) domain="Kinase";
                else if ($3 >= 40351620 && $3 <= 40367741) domain="WD40";
                else domain="UNKNOWN";
                print $6"_UKBB_DP10_"domain, $1"_"$2"_"$4"_"$5;
            }' annotated_coding/LRRK2_all_coding.txt > UKBB_LRRK2.SETID



# # Prepare SETID without G2019S

sed '/40340400/d' UKBB_LRRK2.SETID > UKBB_LRRK2_noG2019S.SETID

# RUN SKAT CAUSE WHY NOT

bash UKBB_run_SKAT_wes.sh 2>&1 | tee 27_OCT_SKAT_LOG.txt



# # NON-SYNONYMOUS VARIANTS
# awk -F "\t" '{ if (NR==1 ||$9 ~ /nonsyn/ ) print $0}'

# # LOF VARIANTS
# awk -F "\t" '{ if (NR==1 || $9 ~ /stop/ || $9 ~ /frame/ || ($9 ~ /intronic_splicing/ && $9 ~ /[+-][1-2]/) ) print $0}'

# # ENCODE VARIANTS
# awk -F "\t" '{ if (NR==1 || $22 != "") print $0}'

# # CADD VARIANTS
# awk -F "\t" '{ if (NR==1 || $15 >= 12.37) print $0}'

