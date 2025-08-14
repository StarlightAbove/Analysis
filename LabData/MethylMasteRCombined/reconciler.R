# A reconciliation script between the normal methylation data and the LMS data, 
# necessary to allow comparative analysis calls from MethylMasteR.

# Copying all normal data and LMS data w/ sheets.
files_to_copy <- list.files("./LabData/Normal_smooth_muscle_EPIC_data/idat_files")
file.copy(from = paste0("./LabData/Normal_smooth_muscle_EPIC_data/idat_files/", files_to_copy), to = "./LabData/MethylMasteRCombined", overwrite = F) 
files_to_copy <- list.files("./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files")
file.copy(from = paste0("./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files/", files_to_copy), to = "./LabData/MethylMasteRCombined", overwrite = F) 

# Combining sheets.
methylSheet <- read.csv("./LabData/MethylMasteRCombined/Sample_Sheet.csv")
normalSheet <- read.csv("./LabData/MethylMasteRCombined/Sample_Sheet_Normal.csv")
SampleSheet <- rbind(methylSheet, normalSheet)
SampleSheet <- SampleSheet  %>% dplyr::select(-c("X"))
write.csv(SampleSheet, file = "./LabData/MethylMasteRCombined/Sample_Sheet_Combined.csv")
