require(EaCoN)
require(readxl)
require(dplyr)
require(readr)
spSheetLMS <- read_excel("LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/design_13LMS_SNP_array.xlsx")
spExtra <- read_excel("LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/design_13LMS_CNVs_other_info_12Aug2025.xlsx")
merged_df <- left_join(spSheetLMS, spExtra, by = "STT")
bth_file <- merged_df %>% dplyr::select("AT Array", "GC Array", "Sample") %>%
  dplyr::rename("ATChannelCel" = "AT Array", "GCChannelCel" = "GC Array", SampleName = Sample) %>%
  dplyr::mutate(ATChannelCel = paste0("~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CEL_files/", ATChannelCel,".CEL"),
                GCChannelCel = paste0("~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CEL_files/", GCChannelCel,".CEL"))
write_tsv(x = bth_file, file = "~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/proc.txt")
setwd("~/Documents/Analysis/LabData/ReprocessedLMS")
ts <- read_tsv("~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/proc.txt")
for(tsx in seq_len(nrow(ts))){
  rw <- ts[tsx, ]
  print(rw)
  OS.Process(ATChannelCel = rw$ATChannelCel, GCChannelCel = rw$GCChannelCel, samplename = rw$SampleName, apt.build = "na33.r4")
}
furtherProcFilesLMS <- list.files(
  path = "~/Documents/Analysis/LabData/ReprocessedLMS",             # Start search from the current working directory
  pattern = "\\.RDS$",    # Regex pattern to match files ending with '.csv'
  recursive = TRUE,       # Search within all subdirectories
  full.names = TRUE       # Return the full relative path to the files
)
for(fpfL in furtherProcFilesLMS){
  Segment.ff(RDS.file = fpfL, segmenter = "ASCAT", BAF.filter = 0.9)
  Segment.ff(RDS.file = fpfL, segmenter = "SEQUENZA", BAF.filter = 0.9)
}

CNVCallsAscat <- list.files(
  path = "~/Documents/Analysis/LabData/ReprocessedLMS",             # Start search from the current working directory
  pattern = "\\.SEG.ASCAT.RDS$",    # Regex pattern to match files ending with '.csv'
  recursive = TRUE,       # Search within all subdirectories
  full.names = TRUE       # Return the full relative path to the files
)
for(CNVCall in CNVCallsAscat){
  ASCN.ff(RDS.file = CNVCall)
}
setwd("~/Documents/Analysis") # On MacBook
