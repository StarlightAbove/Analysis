
sampSheetXLSX <- readxl::read_excel("./LabData/Normal_smooth_muscle_EPIC_data/Uterine_smooth_muscle_5samples_02July2025.xlsx")
# write.csv(sampSheetXLSX, file = "./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files/baseSheet.csv")
fileList <- unique(str_remove(str_remove(string = list.files("./LabData/Normal_smooth_muscle_EPIC_data/idat_files", 
                                                      pattern = "\\.idat$"), pattern = "_Grn.idat"), pattern = "_Red.idat"))
SenPos <- c()
SenID <- c()
for(i in fileList){
  splitPos <- str_split(i, pattern = "_")[[1]]
  SenID <- append(SenID, splitPos[1])
  SenPos <- append(SenPos, splitPos[2])
}
  
IlluminaSheet <- data.frame(Sample_Name = fileList, 
                              Sample_Group = "normal", 
                              Sample_Plate = sampSheetXLSX$Plate, 
                              Pool_ID = NA, 
                              Sentrix_Position = SenPos, 
                              Sentrix_ID = SenID,
                              Basename = fileList, 
                              Platform = "EPIC",
                              Batch = 1)
  
  write.csv(IlluminaSheet, file = "./LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet.csv")
  