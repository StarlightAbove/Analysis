
xlsxSamptoIllumina <- function(dirpath){
  sampSheetXLSX <- readxl::read_excel(dirpath)
  # write.csv(sampSheetXLSX, file = "./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files/baseSheet.csv")
  
  IlluminaSheet <- data.frame(Sample_Name = paste0(sampSheetXLSX$Sentrix_ID, "_", sampSheetXLSX$Sentrix_Position), 
                              Sample_Group = "tumor", 
                              Sample_Plate = sampSheetXLSX$Sample_Plate, 
                              Pool_ID = NA, 
                              Sentrix_Position = sampSheetXLSX$Sentrix_Position, 
                              Sentrix_ID = sampSheetXLSX$Sentrix_ID,
                              Basename = paste0(sampSheetXLSX$Sentrix_ID, "_", sampSheetXLSX$Sentrix_Position), 
                              Platform = "EPIC",
                              Batch = 1)
  
  write.csv(IlluminaSheet, file = "./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files/Sample_Sheet.csv")
  
}

xlsxSamptoIllumina("./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/design_13LMS_EPIC_array.xlsx")
