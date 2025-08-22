# Processing existing .txt to .csv & filtering.
library(tidyverse)

dataList <- paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/",
  list.files(path = "./LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls", pattern = "\\.txt$"))
outputDir <- paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV")

# Simple algorithm that uses the .txt file design to extract the interesting data, and save it as a .csv.
for(dl in dataList){
  fName <- str_split(dl, pattern = "/")[[1]]
  fName <- str_remove(fName[length(fName)], ".txt")
  lines <- readLines(dl)
  lines <- lines[37:(match(T, grepl("MutScore	MutCall", lines)) - 1)]
  writeLines(lines, "trialSubset.txt")
  FASST2Calls <- read.table(file = "trialSubset.txt", sep = "\t")
  names(FASST2Calls) <- FASST2Calls[1, ]
  FASST2Calls <- FASST2Calls %>% dplyr::filter(Event == "CN Loss" | Event == "CN Gain")
  file.remove("trialSubset.txt")
  write.csv(FASST2Calls, paste0(outputDir, "/", fName, ".csv"))
}

# Creation of a matching matrix to associate methylation identifiers to the .csv files.
# Matches CEL files, CSV files, and .idat files in one to prevent extra match processing.
designFileSNP <- readxl::read_excel("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/design_13LMS_SNP_array.xlsx")
designFileSNPExt <- readxl::read_excel(
  "~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/design_13LMS_CNVs_other_info_12Aug2025.xlsx")
designFileEPIC <- readxl::read_excel("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/design_13LMS_EPIC_array.xlsx")
snpfiles <- merge(designFileSNP, designFileSNPExt, by = "STT", all = T)
allfile <- merge(designFileEPIC, snpfiles, by = "STT", all = T) %>% dplyr::select(
  c("STT", "Sample ID", "Sample", "Sentrix_ID", "Sentrix_Position", "AT Array", "GC Array")) %>% 
  dplyr::rename(`CNV_Label` = Sample, Sample_ID = `Sample ID`, GC_Array = `GC Array`, AT_Array = `AT Array`)
write.csv(allfile, paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/correlative.csv"))