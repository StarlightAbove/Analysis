library(GenomicRanges)
library(readxl)

softOutputProcessingLM <- function(outputDir){
  IDATSampleSheet <- read.csv("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
  chasFile <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_ChAS_CNVs.xlsx") %>% select(File, `Mean Log2Ratio`, Chromosome, Type, `Full Location`)
  fasst2file <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_FASST2_CNVs.xlsx")
  
  sesameOutput <- read.csv(outputDir) %>% dplyr::select(c("chrom", "loc.start", 
                                          "loc.end", "seg.mean")) %>% dplyr::mutate(
                                            chrom = str_remove_all(chrom, "chr")) %>% filter(
                                            !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                            chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                            CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                            seg.mean < -0.3 ~ "Deletion",
                                            TRUE ~ "Normal"), type = "SeSAMe") 
  
  caseNameCorrelation <- str_split(outputDir, "_")[[1]]
  sentrixID <- caseNameCorrelation[2]
  sentrixPos <- gsub(".csv", replacement = "", x = caseNameCorrelation[3])
  SUHMatch <- IDATSampleSheet %>% filter(Sentrix_ID == sentrixID, 
                                         Sentrix_Position == sentrixPos) %>% select(SUH) %>% dplyr::pull(SUH)
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
                                                
  
  print(colnames(fasst2filtered))
  print(colnames(chasFiltered))
  print(colnames(sesameOutput))
  rbind(sesameOutput, fasst2filtered, chasFiltered) %>% filter(!is.na(chrom))
}
softOutputProcessingLMS <- function(caseName){
  SampSheet <- read.csv("./cases/SingleFileCase/SampleSheet.csv")
  caseMethyl <- unique(str_remove(str_remove( 
    list.files(paste0("./cases/", caseName, "/maskedMethylation/")), "_Grn.idat"), 
    "_Red.idat"))
  
  outputMethyl <- read.csv(paste0("./Outputs/SeSAMe/LMS/segments_", caseMethyl, ".csv")) %>%
    dplyr::select(c("chrom", "loc.start", 
                    "loc.end", "seg.mean")) %>% dplyr::mutate(
                      chrom = str_remove_all(chrom, "chr")) %>% filter(
                        !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                          chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                            CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                  seg.mean < -0.3 ~ "Deletion",
                                                  TRUE ~ "Normal"), type = "SeSAMe") 
  outputCNV <- read.csv(paste0("./cases/", caseName, "/cnvs/", 
      list.files(paste0("./cases/", caseName, "/cnvs/"), pattern = "\\.csv$"))) %>%
    dplyr::select("Chromosome", "Start", "End", "Segment_Mean") %>%
    dplyr::rename(chrom = "Chromosome", loc.start = "Start", loc.end = "End", seg.mean = "Segment_Mean") %>%
    dplyr::mutate(
      CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                            seg.mean < -0.3 ~ "Deletion",
                            TRUE ~ "Normal"), type = "SNP", chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))
  
  
  
  # Correlating case barcodes to output segments & CNVs.
  rbind(outputMethyl, outputCNV)
}



