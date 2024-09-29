

wd <- getwd()

packrat::init("~/runs/eyu8/library/SKAT")
library(SKAT)
library(MetaSKAT)


setwd(wd)


require(methods)
require(Xmisc)
parser <- ArgumentParser$new()

File.Meta.SKATO = paste0("METASKAT/","METASKAT.robust.results.skato")
File.Meta.BURDEN = paste0("METASKAT/","METASKAT.robust.results.burden")


# Convert the output to a vector of cohort names
all_cohorts <- c("FRENCH","RUSSIA","USA","ISRAEL","UKBB","AMP_PD")

filtered_cohorts <- c()

for (cohort in all_cohorts) {
  File.SetID <- paste0("SKAT/", cohort, ".SETID")

  # Check if File.SetID exists
  if (file.exists(File.SetID)) {
    filtered_cohorts <- c(filtered_cohorts, cohort)
  } else {
    print(paste0(cohort," is excluded as it does not have the SETID file for"))
  }
}


File.Mat.vec<-rep("",6)
File.SetInfo.vec<-rep("",6)

i = 1  
for(cohort in filtered_cohorts){
    File.Bed <- paste0("SKAT/", cohort, ".bed")
      File.Bim <- paste0("SKAT/", cohort, ".bim")
      File.Fam <- paste0("SKAT/", cohort, ".fam")
      File.SetID <- paste0("SKAT/", cohort, ".SETID")
      File.SSD <- paste0("SKAT/", cohort, ".SSD")
      File.Info <- paste0("SKAT/", cohort, ".info")
      File.Mat <- paste0("SKAT/", cohort, ".mat")
      File.SetInfo <- paste0("SKAT/", cohort, ".MInfo")
      File.cov <- paste0("covar_",cohort,".txt")
      File.Results.SKATO <- paste0("SKAT/", cohort, ".results.skato")
      File.Results.BURDEN <- paste0("SKAT/", cohort, ".results.burden")


    Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
    SSD.INFO<-Open_SSD(File.SSD, File.Info)

    FAM<-Read_Plink_FAM_Cov(File.Fam, File.cov, Is.binary =TRUE, cov_header=TRUE)
    y<-FAM$Phenotype
    Age<-FAM$Age
    Sex<-FAM$Sex.y
    N.Sample<-length(y)
    obj<-SKAT_Null_Model(y ~ Sex + Age, out_type="D")

   out.skato<-SKATBinary_Robust.SSD.All(SSD.INFO, obj, method="SKATO")
   out.skato.burden<-SKATBinary_Robust.SSD.All(SSD.INFO, obj, method="Burden")
    write.table(out.skato$results, file=File.Results.SKATO, col.names = TRUE, row.names = FALSE)
    write.table(out.skato.burden$results, file=File.Results.BURDEN, col.names = TRUE, row.names = FALSE)

    re1<-Generate_Meta_Files(obj, File.Bed, File.Bim, File.SetID, File.Mat, File.SetInfo, N.Sample)

    File.Mat.vec[i] <- File.Mat
    File.SetInfo.vec[i] <- File.SetInfo

    i <- i + 1
}


Cohort.Info <- Open_MSSD_File_2Read(File.Mat.vec, File.SetInfo.vec)

out.skato.burden <- MetaSKAT_MSSD_ALL(Cohort.Info, method = "Burden")
out.skato<- MetaSKAT_MSSD_ALL(Cohort.Info, method = "optimal")
write.table(out.skato, file= File.Meta.SKATO, col.names = TRUE, row.names = FALSE)
write.table(out.skato.burden, file= File.Meta.BURDEN, col.names = TRUE, row.names = FALSE)
