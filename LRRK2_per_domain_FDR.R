install.packages("openxlsx")
install.packages("readr")
# Load the required libraries
library(readr)
library(dplyr)
library(tidyr)
library(stats)
library(openxlsx)


dir_path <- "/Users/cemparlar/Desktop/Gan-Or/*JULY2024_yesnoG2019S_LRRK2/RESULTS/noG2019Srerun2/"


# Load your data files
META_LRRK2 <- read.table(paste0(dir_path, "METASKAT.robust.results.skato"), sep = " ", header = TRUE)
UKBB_LRRK2 <- read.table(paste0(dir_path, "UKBB.results.skato"), sep = " ", header = TRUE)
AMP_PD_LRRK2 <- read.table(paste0(dir_path, "AMP_PD.results.skato"), sep = " ", header = TRUE)
McGill_LRRK2 <- read.table(paste0(dir_path, "FRENCH.results.skato"), sep = " ", header = TRUE)
Columbia_LRRK2 <- read.table(paste0(dir_path, "USA.results.skato"), sep = " ", header = TRUE)
Sheba_LRRK2 <- read.table(paste0(dir_path, "ISRAEL.results.skato"), sep = " ", header = TRUE)
Pavlov_LRRK2 <- read.table(paste0(dir_path, "RUSSIA.results.skato"), sep = " ", header = TRUE)


extract_components <- function(data) {
  data <- separate(data, SetID, into = c("Gene", "Domain", "narrow_wide", "Variant_Type"), sep = "~", remove = FALSE)
  return(data)
}

META_LRRK2 <- extract_components(META_LRRK2)
META_LRRK2$p.value <- as.numeric(META_LRRK2$p.value)

UKBB_LRRK2 <- extract_components(UKBB_LRRK2)
UKBB_LRRK2$P.value <- as.numeric(UKBB_LRRK2$P.value)

AMP_PD_LRRK2 <- extract_components(AMP_PD_LRRK2)
AMP_PD_LRRK2$P.value <- as.numeric(AMP_PD_LRRK2$P.value)

McGill_LRRK2 <- extract_components(McGill_LRRK2)
McGill_LRRK2$P.value <- as.numeric(McGill_LRRK2$P.value)

Columbia_LRRK2 <- extract_components(Columbia_LRRK2)
Columbia_LRRK2$P.value <- as.numeric(Columbia_LRRK2$P.value)

Sheba_LRRK2 <- extract_components(Sheba_LRRK2)
Sheba_LRRK2$P.value <- as.numeric(Sheba_LRRK2$P.value)

Pavlov_LRRK2 <- extract_components(Pavlov_LRRK2)
Pavlov_LRRK2$P.value <- as.numeric(Pavlov_LRRK2$P.value)

# Keep only wide for LRRK2

META_LRRK2 <- META_LRRK2 %>% filter(narrow_wide == "wide")
UKBB_LRRK2 <- UKBB_LRRK2 %>% filter(narrow_wide == "wide")
AMP_PD_LRRK2 <- AMP_PD_LRRK2 %>% filter(narrow_wide == "wide")
McGill_LRRK2 <- McGill_LRRK2 %>% filter(narrow_wide == "wide")
Columbia_LRRK2 <- Columbia_LRRK2 %>% filter(narrow_wide == "wide")
Sheba_LRRK2 <- Sheba_LRRK2 %>% filter(narrow_wide == "wide")
Pavlov_LRRK2 <- Pavlov_LRRK2 %>% filter(narrow_wide == "wide")

# Remove narrow_wide column
META_LRRK2 <- subset(META_LRRK2, select = -narrow_wide)
UKBB_LRRK2 <- subset(UKBB_LRRK2, select = -narrow_wide)
AMP_PD_LRRK2 <- subset(AMP_PD_LRRK2, select = -narrow_wide)
McGill_LRRK2 <- subset(McGill_LRRK2, select = -narrow_wide)
Columbia_LRRK2 <- subset(Columbia_LRRK2, select = -narrow_wide)
Sheba_LRRK2 <- subset(Sheba_LRRK2, select = -narrow_wide)
Pavlov_LRRK2 <- subset(Pavlov_LRRK2, select = -narrow_wide)
# Do FDR correction
META_LRRK2$Post_FDR <- p.adjust(META_LRRK2$p.value, method = "BH")
UKBB_LRRK2$Post_FDR <- p.adjust(UKBB_LRRK2$P.value, method = "BH")
AMP_PD_LRRK2$Post_FDR <- p.adjust(AMP_PD_LRRK2$P.value, method = "BH")
McGill_LRRK2$Post_FDR <- p.adjust(McGill_LRRK2$P.value, method = "BH")
Columbia_LRRK2$Post_FDR <- p.adjust(Columbia_LRRK2$P.value, method = "BH")
Sheba_LRRK2$Post_FDR <- p.adjust(Sheba_LRRK2$P.value, method = "BH")
Pavlov_LRRK2$Post_FDR <- p.adjust(Pavlov_LRRK2$P.value, method = "BH")

# Write them into 