library(GenomicRanges)
matchingAlgo <- function(caseName, tech, whichData){
  
  # Pre-processing.
  case <- unique(list.files(path = paste0("./cases/", caseName, "/maskedMethylation"), 
                            pattern = ".idat") %>% str_remove(pattern = "_Red.idat") %>% 
                   str_remove(pattern = "_Grn.idat"))
  actualCaseOutput <- paste0("./Outputs/", tech, "/", whichData, "/", case, ".csv")
  cnv <- list.files(path = paste0("./cases/", caseName, "/cnvs"), pattern = ".csv")
  cnvPath <- paste0("./cases/", caseName, "/cnvs/", cnv)
  
  cnv <- read_csv(cnvPath)
  case <- read_csv(actualCaseOutput)
  
  cnv <- cnv %>% dplyr::select(-c("GDC_Aliquot")) %>%
    mutate(CNVStatus = case_when(cnv$Segment_Mean > 0.3 ~ "Amplification",
                                 cnv$Segment_Mean < -0.3 ~ "Deletion",
                                 TRUE ~ "Normal")) %>% filter(cnv$Num_Probes > 100) %>% 
    dplyr::select(-c("Num_Probes", "...1")) %>% mutate(type = "SNP") %>% dplyr::rename(chrom = Chromosome, loc.start = Start, loc.end = End, seg.mean = Segment_Mean)
  
  # Change labels based on tech + add rename if necessary to match the above 
  # dplyr::rename.
  case <- case %>% dplyr::select(-c("ID", "bstat")) %>% 
    mutate(CNVStatus = case_when(case$seg.mean > 0.2 ~ "Amplification",
                                 case$seg.mean < -0.2 ~ "Deletion",
                                 TRUE ~ "Normal"), 
           chrom = as.numeric(str_remove_all(chrom, "chr"))) %>% 
    dplyr::arrange(chrom) %>% mutate(type = tech) %>%
    dplyr::select(-c("num.mark", "pval", "seg.median", "...1"))
  comb <- rbind(cnv, case)
  comb
}

matchingAlgoLM <- function(outputDir){
  IDATSampleSheet <- read.csv("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
  chasFile <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_ChAS_CNVs.xlsx") %>% select(File, `Mean Log2Ratio`, Chromosome, Type, `Full Location`)
  fasst2file <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_FASST2_CNVs.xlsx")
  
  caseNameCorrelation <- str_split(outputDir, "/")[[1]][[5]]
  caseNameCorrelation <- str_split(caseNameCorrelation, "_")[[1]]
  sentrixID <- caseNameCorrelation[1]
  sentrixPos <- gsub(".csv", replacement = "", x = caseNameCorrelation[2])
  SUHMatch <- IDATSampleSheet %>% filter(Sentrix_ID == sentrixID, 
                                         Sentrix_Position == sentrixPos) %>% select(SUH) %>% dplyr::pull(SUH)
  
  case <- read.csv(outputDir)
  case <- case %>% dplyr::select(-c("ID", "bstat")) %>% 
    mutate(CNVStatus = case_when(case$seg.mean > 0.2 ~ "Amplification",
                                 case$seg.mean < -0.2 ~ "Deletion",
                                 TRUE ~ "Normal"), 
           chrom = as.numeric(str_remove_all(chrom, "chr"))) %>% 
    dplyr::arrange(chrom) %>% mutate(type = "Conumee") %>%
    dplyr::select(-c("num.mark", "pval", "seg.median", "X"))
  
  
  # SNP Filter.
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
  
  rbind(case, chasFiltered, fasst2filtered) %>% dplyr::arrange(chrom) %>% dplyr::filter(!(is.na(chrom)))
}


