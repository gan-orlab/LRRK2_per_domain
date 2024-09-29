#!/bin/bash

module load vcftools
module load StdEnv/2020 plink/1.9b_6.21-x86_64

read BASE_DIR gene cohort_name DP cohort_folder<<< $@
if [[ -z $gene ]]; then echo "ERROR: gene (1st arg) not specified"; exit 42; fi
if [[ -z $cohort_name ]]; then echo "ERROR: cohort_name (2nd arg) not specified"; exit 42; fi
if [[ -z $DP ]]; then echo "ERROR: depth (3rd arg) not specified"; exit 42; fi


ANALYSIS_DIR=$BASE_DIR/$gene/$cohort_name/${DP}x/analysis


vcf=$BASE_DIR/$gene/$cohort_name/${DP}x/$gene.$cohort_name.DP$DP.vcf
sex=$cohort_folder/sex_$cohort_name.txt
pheno=$cohort_folder/pheno_$cohort_name.txt
covar=$cohort_folder/covar_$cohort_name.txt
if [[ ! -s $vcf ]]; then echo "ERROR: vcf file empty or does not exist; name should be $vcf"; exit 42; fi
if [[ ! -s $sex ]]; then echo "ERROR: sex file empty or does not exist; name should be $sex"; exit 42; fi
if [[ ! -s $pheno ]]; then echo "ERROR: pheno file empty or does not exist; name should be $pheno"; exit 42; fi
if [[ ! -s $covar ]]; then echo "ERROR: covar file empty or does not exist; name should be $covar"; exit 42; fi

#name=${gene}_${cohort_name}_DP$DP
name=$BASE_DIR/$gene/$cohort_name/${DP}x/${gene}_${cohort_name}_DP$DP

#make ref txt

awk 'BEGIN{FS=OFS="\t"}{if ($0~/^#/) {next;} print $3,$4}' $vcf > $BASE_DIR/$gene/$cohort_name/${DP}x/REF_ALLELE.txt

#update sex and pheno
plink --vcf $vcf --update-sex $sex --make-bed --allow-no-sex --out $name --output-chr M
plink --bfile $name  --pheno $pheno --make-bed --allow-no-sex --out $name --output-chr M

#filtration steps

plink --bfile $name --missing --allow-no-sex --out ${name}.missingindv --output-chr M
awk 'NR>1 && $5>0.1 {print $2}' ${name}.missingindv.lmiss > ${name}.geno10.snpstoremove
plink --bfile $name --exclude ${name}.geno10.snpstoremove --make-bed --allow-no-sex --out ${name}_geno10 --output-chr M

awk 'NR>1 && $6 >0.1{print $1,$2}' ${name}.missingindv.imiss > ${name}.geno10.indvtoremove
plink --bfile ${name}_geno10 --remove ${name}.geno10.indvtoremove --make-bed --allow-no-sex --out ${name}_geno10_ind10 --output-chr M

plink --bfile ${name}_geno10_ind10 --maf 0.05 --filter-controls --make-bed --out ${name}_geno10_ind10_common
plink --bfile ${name}_geno10_ind10_common --hardy  --allow-no-sex --out ${name}_geno10_ind10_hwe --output-chr M
awk '$9<0.001 {print $2}' ${name}_geno10_ind10_hwe.hwe > ${name}_geno10_ind10_hwe.snpstoremove


plink --bfile ${name}_geno10_ind10 --exclude ${name}_geno10_ind10_hwe.snpstoremove --make-bed --allow-no-sex --out ${name}_geno10_ind10_hwe --output-chr M


plink --bfile ${name}_geno10_ind10_hwe --test-missing --allow-no-sex --out ${name}_geno10_ind10_hwe_testmissing --output-chr M

#Bonferonni cut off 0.05/#of SNPs, see how many get removed
cat ${name}_geno10_ind10_hwe_testmissing.missing | tail -n +2 | awk -v x="$(wc -l ${name}_geno10_ind10_hwe.bim)" '{ if ($5<0.05/x) {print $2}}'  > ${name}_geno10_ind10_hwe_testmissing.snpstoremove
#awk 'NR>1 && $5<0.0125 {print $2}' ${name}_geno10_ind10_hwe_testmissing.missing > ${name}_geno10_ind10_hwe_testmissing.snpstoremove

plink --bfile ${name}_geno10_ind10_hwe  --exclude ${name}_geno10_ind10_hwe_testmissing.snpstoremove --make-bed --allow-no-sex --out ${name}_geno10_ind10_hwe_testmiss --output-chr M

covar_name="Sex,Age"

if [ $cohort_name == newyork ]
then
    covar_name="Sex,Age,Ethn"
fi

plink --bfile ${name}_geno10_ind10_hwe_testmiss --a2-allele $BASE_DIR/$gene/$cohort_name/${DP}x/REF_ALLELE.txt --logistic hide-covar --covar $covar --covar-name Age, Sex --ci 0.95 --allow-no-sex  --out $name --output-chr M

plink --bfile ${name}_geno10_ind10_hwe_testmiss --a2-allele $BASE_DIR/$gene/$cohort_name/${DP}x/REF_ALLELE.txt --assoc fisher --out $name --output-chr M

cp ${name}_geno10_ind10_hwe_testmiss* $ANALYSIS_DIR/

cut -d ' ' -f1 ${name}_geno10_ind10_hwe_testmiss.fam > ${name}.final.samples

plink --bfile ${name}_geno10_ind10_hwe_testmiss --a2-allele $BASE_DIR/$gene/$cohort_name/${DP}x/REF_ALLELE.txt --recode vcf --out ${name}_geno10_ind10_hwe_testmiss --output-chr M

paste ${name}.assoc.logistic ${name}.assoc.fisher | sed 's/ \+/\t/g' | awk '{if (NR>1){tmp=$4;$4=$19;$19=tmp} print $1,$2,$3,$4,$19,$6,$7,$9,$10,$12,$17,$18}' OFS="\t" > $ANALYSIS_DIR/$gene.$cohort_name.DP$DP.stats

