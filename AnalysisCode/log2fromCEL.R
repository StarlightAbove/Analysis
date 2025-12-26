# To analyze log2 ratio from CEL files, to provide dataset blanks in existing 
# dataset.
library(affy)
library(oligo)
# Requires copynumber, sequenza
library(oncoscanR)
setwd("/home/Eliza/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CEL_files")
celFiles <- list.celfiles(full.names = TRUE)
rawData <- read.celfiles(celFiles)
print(exprs_data)
setwd("/home/Eliza/Work/Analysis")