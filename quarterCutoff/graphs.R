

# Add anonymized labels for cases
correlative <- read.csv("~/projects/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv")
tags <- paste0("LMS", 1:length(correlative$X))
correlative <- correlative %>% dplyr::mutate(LMSTags = tags)
write.csv(x = correlative, file = "~/projects/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv")
