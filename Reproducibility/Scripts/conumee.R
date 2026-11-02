library(tidyverse)
library(R.utils)
library(conumee2)
library(reshape2)

# Bin sizes (bp) evaluated throughout the analysis.
binLength <- c(10000, 50000, 1e+05, 1e+06)

con2 <- function(locationOfTest, locationOfControl, arrayType, binSize, OutputFolder){
  RGset <- read.metharray.exp(locationOfTest, recursive = T, verbose = T)
  MSet <- preprocessIllumina(RGset)
  CSet <- read.metharray.exp(locationOfControl, recursive = T, verbose = T)
  CSet <- preprocessIllumina(CSet)
  anno <- CNV.create_anno(array_type = c("450k", "EPIC"), bin_minsize = binSize)
  load.data <- CNV.load(MSet)
  load.controls <- CNV.load(CSet)
  x <- CNV.fit(load.data, ref = load.controls, anno)
  x <- CNV.bin(x)
  # x <- CNV.detail(x) (Can be added if there is certain focus on certain loci)
  x <- CNV.segment(x)
  segments <- CNV.write(x, what = "segments")
  bins <- CNV.write(x, what = "bins")
  
  dir.create(paste0("./Outputs/Conumee/",OutputFolder,"/",binSize), recursive = TRUE, showWarnings = FALSE)
  write.csv(bins, paste0("./Outputs/Conumee/",OutputFolder,"/",binSize,"/bins.csv"))
  for(segment in segments){
    write.csv(segment, paste0("./Outputs/Conumee/",OutputFolder,"/",binSize,"/",segment[["ID"]][1], ".csv"))
  }
}

for(bin in binLength){
  con2(locationOfTest = "./LabData/LMS_SNP_EPIC_array_data", 
       locationOfControl = "./LabData/Normal_smooth_muscle_EPIC_data", 
       arrayType = "EPIC",
       binSize = bin,
       OutputFolder = "LabLMS")
} # LabLMS
for(bin in binLength){
  con2(locationOfTest = "./LabData/Normal_smooth_muscle_EPIC_data", 
       locationOfControl = "./Controls/EPICControls/", 
       arrayType = "EPIC",
       binSize = bin,
       OutputFolder = "LabNormals")
} # LabNormals. GEO adjacent-myometrium controls (GSE192354) are used here only
  # to quality-check the lab normals, which are the reference for both tumour cohorts.
for(bin in binLength){
  con2(locationOfTest = "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files", 
       locationOfControl = "./LabData/Normal_smooth_muscle_EPIC_data", 
       arrayType = "EPIC",
       binSize = bin,
       OutputFolder = "LMData/bins") # path read by LMStt() in Analytic Code
} # LM.