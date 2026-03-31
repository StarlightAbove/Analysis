library(sesame)
library(sesameData)
library(reshape2)
idat_dir <- "./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files"
sdfs <- openSesame(idat_dir, func = NULL)

i <- 1
for(sdf in sdfs){
  segs <- sdf
  segments <- cnSegmentation(segs)
  segmentalSignals <- segments[["seg.signals"]]
  bins <- segments[["bin.signals"]]
  bins <- reshape2::melt(bins)
  write.csv(bins, paste0("./Outputs/SeSAMe/LabLMS/bins_",names(sdfs[i]),".csv"), row.names = TRUE)
  write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/LabLMS/segments_", 
                                     names(sdfs[i]), ".csv"), row.names = FALSE)
  i <- i + 1
}

