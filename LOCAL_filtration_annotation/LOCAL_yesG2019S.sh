#!/bin/bash


CrossMap vcf ~/runs/emsom/softwares/liftOver/hg19ToHg38.over.chain.gz FRENCH.DP30.all.genes.vcf hg38.fa FRENCH_hg38.vcf
CrossMap vcf ~/runs/emsom/softwares/liftOver/hg19ToHg38.over.chain.gz ISRAEL.DP30.all.genes.vcf hg38.fa ISRAEL_hg38.vcf
CrossMap vcf ~/runs/emsom/softwares/liftOver/hg19ToHg38.over.chain.gz USA.DP30.all.genes.vcf hg38.fa USA_hg38.vcf
CrossMap vcf ~/runs/emsom/softwares/liftOver/hg19ToHg38.over.chain.gz RUSSIA.DP30.all.genes.vcf hg38.fa RUSSIA_hg38.vcf

plink --vcf "FRENCH_hg38.vcf" --vcf-half-call m --make-bed --out "FRENCH_all_variants" --output-chr M
plink --vcf "ISRAEL_hg38.vcf" --vcf-half-call m --make-bed --out "ISRAEL_all_variants" --output-chr M
plink --vcf "USA_hg38.vcf" --vcf-half-call m --make-bed --out "USA_all_variants" --output-chr M
plink --vcf "RUSSIA_hg38.vcf" --vcf-half-call m --make-bed --out "RUSSIA_all_variants" --output-chr M

for cohort in "FRENCH" "USA" "RUSSIA" "ISRAEL"; do
	plink --bfile "${cohort}_all_variants"  --update-sex cohort_${cohort}/sex_${cohort}.txt --pheno cohort_${cohort}/pheno_${cohort}.txt --pheno-name Status --make-bed --allow-no-sex --keep cohort_${cohort}/pheno_${cohort}.txt --out ${cohort}_formatted_all_variants --output-chr M
done

plink --bfile "FRENCH_formatted_all_variants" --max-maf 0.01 --write-snplist --out FRENCH_rare_variants
plink --bfile "ISRAEL_formatted_all_variants" --max-maf 0.01 --write-snplist --out ISRAEL_rare_variants
plink --bfile "USA_formatted_all_variants" --max-maf 0.01 --write-snplist --out USA_rare_variants
plink --bfile "RUSSIA_formatted_all_variants" --max-maf 0.01 --write-snplist --out RUSSIA_rare_variants

# echo "rs34637584:G:A" >> FRENCH_rare_variants.snplist
# echo "rs34637584:G:A" >> RUSSIA_rare_variants.snplist
echo "rs34637584:G:A" >> ISRAEL_rare_variants.snplist
echo "rs34637584:G:A" >> USA_rare_variants.snplist

# for cohort in "FRENCH" "USA" "RUSSIA" "ISRAEL"; do
# grep "rs34637584" ${cohort}_formatted_all_variants.bim
# done

for cohort in "FRENCH" "USA" "RUSSIA" "ISRAEL"; do
	plink --bfile ${cohort}_formatted_all_variants --extract ${cohort}_rare_variants.snplist --make-bed --out ${cohort}
done

mkdir -p annotation

for cohort in "FRENCH" "USA" "RUSSIA" "ISRAEL"; do
	bgzip -c ${cohort}_hg38.vcf > ${cohort}_hg38.vcf.gz
	tabix -p vcf ${cohort}_hg38.vcf.gz
	awk '{print $1"\t"$4-1"\t"$4"\t"$2}' ${cohort}.bim > ${cohort}_variants_list.bed
	bcftools view -R ${cohort}_variants_list.bed ${cohort}_hg38.vcf.gz -o ${cohort}_hg38_filtered.vcf
	perl ~/runs/sitkicem/annovar/convert2annovar.pl --format vcf4 "${cohort}_hg38_filtered.vcf" --allsample --withfreq --outfile "annotation/${cohort}_recode_convert"
	perl ~/runs/sitkicem/annovar/table_annovar.pl "annotation/${cohort}_recode_convert" ~/runs/sitkicem/annovar/humandb/ --buildver hg38 --out "annotation/${cohort}_withgnomad_clinvar.annovar" --remove --protocol refGene,ljb26_all,dbnsfp47c,clinvar_20240611,gnomad41_genome,avsnp151 --operation g,f,f,f,f,f --nastring .
done

for cohort in "FRENCH" "USA" "RUSSIA" "ISRAEL"; do
perl ~/runs/sitkicem/annovar/table_annovar.pl "annotation/${cohort}_recode_convert" ~/runs/sitkicem/annovar/humandb/ --buildver hg38 --out "annotation/${cohort}_withgnomad_clinvar.annovar" --remove --protocol refGene,ljb26_all,dbnsfp47c,clinvar_20240611,gnomad41_genome,avsnp151 --operation g,f,f,f,f,f --nastring .
done


# perl ~/runs/sitkicem/annovar/table_annovar.pl "annotation/UKBB_recode_convert" ~/runs/sitkicem/annovar/humandb/ --buildver hg38 --out "annotation/UKBB_withgnomad_clinvar.annovar" --remove --protocol refGene,ljb26_all,dbnsfp42c,clinvar_20240611,gnomad41_genome,avsnp151 --operation g,f,f,f,f,f --nastring .
# perl ~/runs/sitkicem/annovar/table_annovar.pl "annotation/AMP_PD.avinput" ~/runs/sitkicem/annovar/humandb/ --buildver hg38 --out "annotation/AMP_PD_withgnomad_clinvar.annovar" --remove --protocol refGene,ljb26_all,dbnsfp42c,clinvar_20240611,gnomad41_genome,avsnp151 --operation g,f,f,f,f,f --nastring .

for cohort in "FRENCH" "ISRAEL" "USA" "RUSSIA"; do
	# NON-SYNONYMOUS VARIANTS
	awk -F "\t" '{ if (NR==1 ||$9 ~ /nonsyn/ ) print $1"_"$2"_"$4"_"$5}' annotation/${cohort}_recode_convert.annovar.hg38_multianno.txt > annotation/${cohort}_nonsyn.txt

	# LOF VARIANTS
	awk -F "\t" '{ if (NR==1 || $9 ~ /stop/ || $9 ~ /frame/ || ($9 ~ /intronic_splicing/ && $9 ~ /[+-][1-2]/) ) print $1"_"$2"_"$4"_"$5}' annotation/${cohort}_recode_convert.annovar.hg38_multianno.txt > annotation/${cohort}_LOF.txt

	# ENCODE VARIANTS
	awk -F "\t" '{ if (NR==1 || $6 == "exonic") print $1"_"$2"_"$4"_"$5}' annotation/${cohort}_recode_convert.annovar.hg38_multianno.txt > annotation/${cohort}_ENCODE.txt

	# CADD VARIANTS
	awk -F "\t" '{ if (NR==1 || $31 >= 20) print $1"_"$2"_"$4"_"$5}' annotation/${cohort}_recode_convert.annovar.hg38_multianno.txt > annotation/${cohort}_CADD.txt
done


## FOR TABLE PREP, PREPARE 
