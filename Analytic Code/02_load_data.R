### ============================================================
### 2. READ ALL DATA
### ============================================================
# Case ID lists, bin sizes, and hg19 reference sizes.
# Requires 00_setup.R (read_csv, getChromInfoFromUCSC).

# Case ID lists for each cohort, the bin sizes analyzed throughout, and
# the hg19 reference sizes used by several of the metric functions above.

# Getting the list of cases
LMS_correlative <- read_csv(paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/correlative.csv"))
LMS_cases <- LMS_correlative$STT
LMS_tags <- setNames(LMS_correlative$LMSTags, LMS_correlative$STT)

LM_SampSheet <- read_csv(paste0(getwd(),"/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv"))
LM_cases <- LM_SampSheet$STT
LM_tags <- setNames(LM_SampSheet$LMSTags, LM_SampSheet$STT)

Normal_SampleSheet <- read_csv("LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")
Normals <- Normal_SampleSheet$Basename
Normal_tags <- setNames(Normal_SampleSheet$NormTags, Normal_SampleSheet$Basename)

bins <- c(10000, 50000, 1e+05, 1e+06)

#### 2a. Prelim Data ----
hg19_info <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(assembled == "TRUE")
hg19_info <- hg19_info[1:22,]
hg19_total <- sum(hg19_info$size)

