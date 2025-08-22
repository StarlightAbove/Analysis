# LMS.
lmsProcessing <- function(caseName){
  sampSheet <- read.csv("./cases/SingleFileCase/SampleSheet.csv") %>% select("Basename", "caseNames")
  match <- sampSheet %>% filter(caseNames == caseName) %>% dplyr::pull("Basename")
  
  # Processing methylation data.
  cnvMethyl <- read.csv(paste0("./Outputs/MethylMaster/LMS/", match, "/autocorrected_regions.csv")) %>%
    dplyr::select(c("Chromosome", "bp.Start", "bp.End", "Mean")) %>% dplyr::rename(
      chrom = "Chromosome",
      loc.start = "bp.Start",
      loc.end = "bp.End",
      seg.mean = "Mean"
    ) %>% dplyr::mutate(
      CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                            seg.mean < -0.3 ~ "Deletion",
                            TRUE ~ "Normal"), type = "MethylMaster"
    )
  cnvMethyl$chrom <- as.numeric(str_remove_all(cnvMethyl$chrom, pattern = "chr"))
  cnvMethyl <- cnvMethyl %>% filter(!(is.na(chrom))) %>% arrange(chrom)
  
  
  
  
  cnvSNP <- read.csv(paste0("./cases/", caseName, "/cnvs/", 
                            list.files(paste0("./cases/", caseName, "/cnvs/"), 
                                       pattern = "\\.csv$"))) %>% dplyr::select("Chromosome", "Start", "End", "Segment_Mean") %>%
    dplyr::rename(chrom = "Chromosome", loc.start = "Start", loc.end = "End", seg.mean = "Segment_Mean") %>%
    dplyr::mutate(
      CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                            seg.mean < -0.3 ~ "Deletion",
                            TRUE ~ "Normal"), type = "SNP", chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))
  
  rbind(cnvMethyl, cnvSNP)
}

# LM.
lmProcessing <- function(outputDir){
  IDATSampleSheet <- read.csv("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
  chasFile <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_ChAS_CNVs.xlsx") %>% select(File, `Mean Log2Ratio`, Chromosome, Type, `Full Location`)
  fasst2file <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_FASST2_CNVs.xlsx")
  
  caseNameCorrelation <- str_split(str_split(outputDir, "/")[[1]][5], pattern = "_")[[1]]
  sentrixID <- caseNameCorrelation[[1]]
  sentrixPos <- caseNameCorrelation[[2]]
  SUHMatch <- IDATSampleSheet %>% filter(Sentrix_ID == sentrixID, 
                                        Sentrix_Position == sentrixPos) %>% select(SUH) %>% dplyr::pull(SUH)
  cnvMethyl <- read.csv(outputDir) %>%
    dplyr::select(c("Chromosome", "bp.Start", "bp.End", "Mean")) %>% dplyr::rename(
      chrom = "Chromosome",
      loc.start = "bp.Start",
      loc.end = "bp.End",
      seg.mean = "Mean"
    ) %>% dplyr::mutate(
      CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                            seg.mean < -0.3 ~ "Deletion",
                            TRUE ~ "Normal"), type = "MethylMaster"
    )
  cnvMethyl$chrom <- as.numeric(str_remove_all(cnvMethyl$chrom, pattern = "chr"))
  cnvMethyl <- cnvMethyl %>% filter(!(is.na(chrom))) %>% arrange(chrom)
  
  # SNP Processing
  chasFiltered <- chasFile[str_detect(chasFile$File, SUHMatch), ]
  chasFiltered <-  tidyr::separate_wider_delim(tidyr::separate_wider_delim(chasFiltered, `Full Location`, names = c("chromosome", "loc"), delim = ":"), `loc`,
                                               names = c("loc.start", "loc.end"), delim = "-") %>% select(-c("chromosome")) %>% dplyr::rename(CNVStatus = "Type", seg.mean = "Mean Log2Ratio", chrom = "Chromosome") %>%
    mutate(
      type = "SNP", 
      CNVStatus = case_when(CNVStatus == "Loss" ~ "Deletion", 
                            CNVStatus == "Gain" ~ "Amplification", 
                            TRUE ~ "Normal")
    ) %>% select(-c("File"))
  fasst2filtered <- fasst2file[str_detect(fasst2file$Sample, SUHMatch), ] %>%
    select(c("Probe Median", "Event", "Chromosome Region"))
  fasst2filtered <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(fasst2filtered, `Chromosome Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
                                                names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") %>% dplyr::rename(CNVStatus = "Event") %>%
    mutate(CNVStatus = case_when(
      CNVStatus == "CN Loss" ~ "Deletion", 
      CNVStatus == "CN Gain" ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% dplyr::rename(seg.mean = "Probe Median")
  fasst2filtered$loc.start <- as.numeric(str_replace_all(fasst2filtered$loc.start, ",", ""))
  fasst2filtered$loc.end <- as.numeric(str_replace_all(fasst2filtered$loc.end, ",", ""))
  fasst2filtered$chrom <- as.numeric(str_replace_all(fasst2filtered$chrom, "chr", ""))
  chasFiltered$chrom <- as.numeric(chasFiltered$chrom)
  chasFiltered$loc.start <- as.numeric(chasFiltered$loc.start)
  chasFiltered$loc.end <- as.numeric(chasFiltered$loc.end)
  chasFiltered$seg.mean <- as.numeric(chasFiltered$seg.mean)
  
  rbind(chasFiltered, fasst2filtered, cnvMethyl) %>% arrange(chrom) %>% dplyr::filter(!(is.na(chrom)))
}

labLMSProc <- function(outputDir){
  # Correlate by STT information between methylation Sentrix and SNP data.
  
}



