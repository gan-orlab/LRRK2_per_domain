#!/bin/bash

module load nixpkgs/16.09
module load gcc/7.3 r/3.5.2

plink --bfile UKBB  --update-sex ~/runs/sitkicem/NEW_ALL_GENES/covar_UKBB_no_proxy/sex_UKBB.txt --pheno ~/runs/sitkicem/NEW_ALL_GENES/covar_UKBB_no_proxy/pheno_UKBB.txt --pheno-name Status --make-bed --allow-no-sex --keep ~/runs/sitkicem/NEW_ALL_GENES/covar_UKBB_no_proxy/pheno_UKBB.txt --out UKBB --output-chr M

srun --mem=40G --cpus-per-task=10 --time=1:00:00 Rscript MetaSKAT_LRRK2_RERUN.R 2>&1|tee METASKAT_LRRK2_RERUN.log
