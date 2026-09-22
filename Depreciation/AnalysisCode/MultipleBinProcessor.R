library(tidyverse)
library(R.utils)
library(conumee2)
library(sesame)
library(sesameData)
library(reshape2)
binLength <-  c(10000,25000,50000,75000,100000, 5e+05, 1e+06)

# Conumee analysis across given bin lengths.
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
} # LabNormals.
for(bin in binLength){
  con2(locationOfTest = "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files", 
       locationOfControl = "./Controls/EPICControls/", 
       arrayType = "EPIC",
       binSize = bin,
       OutputFolder = "LMData")
} # LM.


sdfs <- openSesame("./LabData/LMS_SNP_EPIC_array_data", func = NULL)
for(bin in binLength){
  dir.create(paste0("./Outputs/SeSAMe/LabLMS/bins/",bin))
  dir.create(paste0("./Outputs/SeSAMe/LabLMS/",bin))
  i <- 1
  for(sdf in sdfs){
    segs <- sdf
    segments <- cnSegmentation(segs, tilewidth = bin)
    segmentalSignals <- segments[["seg.signals"]]
    bins <- segments[["bin.signals"]]
    bins <- reshape2::melt(bins)
    write.csv(bins, paste0("./Outputs/SeSAMe/LabLMS/bins/",bin,"/bins_",
                           names(sdfs[i]),".csv"), row.names = TRUE)
    write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/LabLMS/",bin,"/segments_", 
                                       names(sdfs[i]), ".csv"), row.names = FALSE)
    i <- i + 1
  }
} # LabLMS.

sdfs <- openSesame("./LabData/Normal_smooth_muscle_EPIC_data", func = NULL)
for(bin in binLength){
  i <- 1
  for(sdf in sdfs){
    segments <- cnSegmentation(sdf, tilewidth = bin)
    segmentalSignals <- segments[["seg.signals"]]
    bins <- segments[["bin.signals"]]
    bins <- reshape2::melt(bins)
    write.csv(bins, paste0("./Outputs/SeSAMe/Normals/bins/",bin,"/bins_",
                           names(sdfs[i]),".csv"), row.names = TRUE)
    write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/Normals/",bin,"/segments_", 
                                       names(sdfs[i]), ".csv"), row.names = FALSE)
    print(i)
    i <- i + 1
  }
} # Lab Normals.
rm(bins, sdf, segmentalSignals, segments, bin, segs)

sdfs <- openSesame("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files", func = NULL)
binLength <- c(100000, 5e+05, 1e+06)
for(bin in binLength){
  i <- 1
  for(sdf in sdfs){
    segments <- cnSegmentation(sdf, tilewidth = bin)
    segmentalSignals <- segments[["seg.signals"]]
    bins <- segments[["bin.signals"]]
    bins <- reshape2::melt(bins)
    write.csv(bins, paste0("./Outputs/SeSAMe/LM/bins/",bin,"/bins_",
                           names(sdfs[i]),".csv"), row.names = TRUE)
    write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/LM/bins/",bin,"/segments_", 
                                       names(sdfs[i]), ".csv"), row.names = FALSE)
    print(i)
    i <- i + 1
  }
} # LM.
rm(bins, sdf, segmentalSignals, segments, bin, segs)