

# Add anonymized labels for cases
correlative <- read.csv("~/projects/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv")
tags <- paste0("LMS", 1:length(correlative$X))
correlative <- correlative %>% dplyr::mutate(LMSTags = tags)
write.csv(x = correlative, file = "~/projects/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv")

correlative <- read.csv("~/projects/Analysis/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
tags <- paste0("LM", 1:length(correlative$X))
correlative <- correlative %>% dplyr::mutate(LMTags = tags)
write.csv(x = correlative, file = "~/projects/Analysis/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")

Sample_Sheet_Normal <- read_csv("~/projects/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")
tags <- paste0("NORM", 1:length(Sample_Sheet_Normal$...1))
Sample_Sheet_Normal <- Sample_Sheet_Normal %>% dplyr::mutate(NormTags = tags)
write.csv(x = Sample_Sheet_Normal, file = "~/projects/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")
