library("TCGAbiolinks")
library("rjson")
library("readr")

# Selected Cohort Characteristics:
# Demographic: White Non-Hispanic Females
# Project: TCGA-SARC OPEN
# Data Types: Analytical on DNA Methylation (Illumina EPIC), Control on Preprocessed CNVs (Affymetrix SNP 6.0)
# Cancer: Leiomyosarcoma NOS
# Location: Retroperitoneum & Peritoneum
# Vitals: Alive
# 10 sample patients were selected as available when filters were applied. 

jsonq <- rjson::fromJSON(file ="./BaseData/cohort.2025-06-12.json") # Cohort IDs. Change json file to json downloaded from the manifest.
barcodes <- NULL;
for(p in jsonq){
  tmp <- p$submitter_id
  barcodes <- rbind(barcodes, tmp)
}
query <- TCGAbiolinks::GDCquery(project = "TCGA-SARC", access = "open", data.category = c("DNA Methylation"), barcode = barcodes, data.type = "Methylation Beta Value")
queryCNV <- TCGAbiolinks::GDCquery(project = "TCGA-SARC", access = "open", data.category = c("Copy Number Variation"), barcode = barcodes, sample.type = "Primary Tumor", data.type = "Copy Number Segment")
MaskedMethyl <- TCGAbiolinks::GDCquery(project = "TCGA-SARC", access = "open", data.category = c("DNA Methylation"), barcode = barcodes, sample.type = "Primary Tumor", data.type = "Masked Intensities")
TCGAbiolinks::GDCdownload(query = query, method = "api", directory = "methylQuery")
TCGAbiolinks::GDCdownload(query = queryCNV, method = "api", directory = "CNVQuery")
TCGAbiolinks::GDCdownload(query = MaskedMethyl, method = "api", directory = "MaskedMethyl")


# Collation.
# Binning and collating each dataset to a single case.
caseNums <- unique(MaskedMethyl[[1]][[1]]$cases.submitter_id)
dir.create("cases")
for(csns in caseNums){
  dir.create(file.path("cases", csns, "cnvs"), recursive = TRUE)
  dir.create(file.path("cases", csns, "methylation"), recursive = TRUE)
  dir.create(file.path("cases", csns, "maskedMethylation"), recursive = TRUE)
}

# Cleaning environment.
rm(csns)
rm(tmp)
rm(p)
rm(jsonq)
rm(barcodes)

# Sorting all files into correct correlated cases for easy comparison.
bspthcnvq <- "./BaseData/CNVQuery/TCGA-SARC/Copy_Number_Variation/Copy_Number_Segment/"
bspthmethyl <- "./BaseData/methylQuery/TCGA-SARC/DNA_Methylation/Methylation_Beta_Value/"
bspthmasked <- "./BaseData/MaskedMethyl/TCGA-SARC/DNA_Methylation/Masked_Intensities/"
bspthcases <- "./cases/"

# CNVs and Methylation into correct cases.
for(i in caseNums){
  colIDs <- queryCNV[[1]][[1]][which(queryCNV[[1]][[1]]$cases.submitter_id == i), ]$id
  movePath <- paste0(bspthcnvq, colIDs)
  toGoPath <- paste0(bspthcases, i, "/cnvs/")
  file.copy(movePath, toGoPath, recursive = TRUE, overwrite = TRUE)
}
for(i in caseNums){
  colIDs <- query[[1]][[1]][which(query[[1]][[1]]$cases.submitter_id == i), ]$id
  movePath <- paste0(bspthmethyl, colIDs)
  toGoPath <- paste0(bspthcases, i, "/methylation/")
  file.copy(movePath, toGoPath, recursive = TRUE)
}
for(i in caseNums){
  colIDs <- MaskedMethyl[[1]][[1]][which(MaskedMethyl[[1]][[1]]$cases.submitter_id == i), ]$id
  movePath <- paste0(bspthmasked, colIDs)
  toGoPath <- paste0(bspthcases, i, "/maskedMethylation/")
  file.copy(movePath, toGoPath, recursive = TRUE, overwrite = TRUE)
}

# Pre-processing identifiers done.
