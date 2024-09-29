
#salloc -c 40 --mem=160G 

library(data.table)

#Read table with all fields
ukb <- as.data.frame(fread("/project/rpp-aevans-ab/neurohub/ukbb/tabular/current.csv"))

#Write field of interest
field <- c("20002","41270","20111","20110","20107","22001","22006","22009","22000","189","34","21022","22021","22019","22027")

#Change into pattern recognisable by grep
pattern <- paste0("^",field,"-",collapse = "|")

#Selct field of interest
ukb_filtered <- ukb[,c(1,grep(pattern,names(ukb)))]

#Select PD self-reported (1262) and ICD10 code (32330):
#SELF-REPORT DATA
ukb_self <- ukb_filtered[,c(1,grep("^20002-",names(ukb_filtered)))]

PD_self <- ukb_self[which(apply(ukb_self[,-1],1,function(i){any(i == 1262)})),]$eid

#ICD10 DATA
ukb_ICD10 <- ukb_filtered[,c(1,grep("^41270-",names(ukb_filtered)))]


PD_ICD <- ukb_ICD10[which(apply(ukb_ICD10[,-1],1,function(i){any(i == "G20")})),]$eid



#Example Code: PD <- union(PD_self, PD_ICD)
PD <- union(PD_self,PD_ICD)

#Select proxy cases from father, mother, sibling, and exclude PD cases
ukb_proxy <- ukb_filtered[,c(1,grep("^20111-|^20110-|^20107-",names(ukb_filtered)))]
PD_proxy_notexcluding_PD <- ukb_proxy[which(apply(ukb_proxy,1,function(i){any(i == 11)})),]$eid
PD_proxy <- setdiff(PD_proxy_notexcluding_PD, PD)

#Select the rest as controls
PD_control <- setdiff(ukb_filtered$eid,union(PD_proxy,PD))

#Perform filter for samples with known issue (aneupleudy, missingness, het outlier) and relatedness (0 = no closer than 3rd degree relative) & ancestry filter (1 = causacian)
ukb_unrelated <- readLines("~/runs/go_lab/GRCh37/ukbb/ukbb_raw_data_no_cousins.txt")

ukb_filtered_unrelated <- ukb_filtered[ukb_filtered$eid %in% ukb_unrelated,]
ukb_filtered_unrelated_euro <- ukb_filtered[ukb_filtered$"22006-0.0" %in% 1,]
ukb_filtered_unrelated_euro_aneu <- ukb_filtered_unrelated_euro[!(ukb_filtered_unrelated_euro$"22019-0.0" %in% 1),]
ukb_filtered_unrelated_euro_aneu_miss <- ukb_filtered_unrelated_euro_aneu[!(ukb_filtered_unrelated_euro_aneu$"22027-0.0" %in% 1),]

#Add pc's
PC <- as.data.frame(fread("~/runs/go_lab/GRCh37/ukbb/pc_euro.txt"))
PC$IID <- NULL
names(PC)[1] <- "eid"

#Select covariates used in GWAS
covar_field <- c("22001","22009","22000","189","34","21022")
covar <- paste0("^",covar_field,"-",collapse = "|")
ukb_covar <- ukb_filtered_unrelated_euro_aneu_miss[,c(1,grep(covar,names(ukb_filtered_unrelated_euro_aneu_miss)))]
ukb_covar_pc10 <- merge(ukb_covar[,1:6],PC)
#Rename covariates
names(ukb_covar_pc10) <- c("ID", "YearAtBirth", "Townsend", "AgeAtRecruit", "Batch", "Sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10")

#Save seperate covariate for each type of disorder(Tic, OCD, ADHD) and their respective  case and control

write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% PD,], "~/runs/sitkicem/OCT_UKBB/ukbb_PD_covar.txt", quote = F, row.names = F)
write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% PD_proxy,], "~/runs/sitkicem/OCT_UKBB/ukbb_PD_proxy_covar.txt", quote = F, row.names = F)
write.csv(ukb_covar_pc10[ukb_covar_pc10$ID %in% PD_control,], "~/runs/sitkicem/OCT_UKBB/ukbb_PD_control_covar.txt", quote = F, row.names = F)
