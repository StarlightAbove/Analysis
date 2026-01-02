labNmrlProc <- function(Sentrix, Technology, binSize){
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(
    "~/Work/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv") %>% filter(Basename == Sentrix)
  
  # cnvMatch <- read.csv(paste0("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV/", 
  # correlationSheet$CNV_Label[1], "_events.csv"))
  
  # methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
  #paste0(correlationSheet$Sentrix_ID[1],"_", 
  #     correlationSheet$Sentrix_Position[1], "/"), 
  # "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0("~/Work/Analysis/Outputs/Conumee/LabNormals/", binSize,"/", paste0(Sentrix,".csv")))
  
  sesameMatch <- read.csv(paste0("~/Work/Analysis/Outputs/SeSAMe/Normals/", binSize,"/", paste0("segments_",Sentrix,".csv")))
  
  # cnvMatch <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(cnvMatch, `Chromosome.Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
  #                                        names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") %>% dplyr::rename(CNVStatus = "Event") %>%
  # mutate(CNVStatus = case_when(
  #   CNVStatus == "CN Loss" ~ "Deletion", 
  #   CNVStatus == "CN Gain" ~ "Amplification", 
  #   TRUE ~ "Normal"
  # )) %>% dplyr::rename(seg.mean = "Probe.Median") %>% 
  #  dplyr::select(c(chrom, loc.start, loc.end, CNVStatus, seg.mean, type)) %>% mutate(loc.start = as.numeric(gsub(",","",loc.start)), loc.end = as.numeric(gsub(",","",loc.end)))
  #cnvMatch$chrom <- as.numeric(str_replace_all(cnvMatch$chrom, "chr", ""))
  #cnvMatch <- cnvMatch %>% filter(!is.na(chrom))
  
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.2 ~ "Deletion", 
        seg.mean >= 0.2 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- methylMatch %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- conumeeMatch %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                 seg.mean < -0.3 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- sesameMatch %>% arrange(chrom)
  }
  
  
  
  combinedSet
}
lmProcessingSes <- function(outputDir){
  IDATSampleSheet <- read.csv("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
  chasFile <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_ChAS_CNVs.xlsx") %>% dplyr::select(File, `Mean Log2Ratio`, Chromosome, Type, `Full Location`)
  fasst2file <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_FASST2_CNVs.xlsx")
  
  caseNameCorrelation <- str_split(str_split(outputDir, "/")[[1]][7], pattern = "_")[[1]]
  sentrixID <- caseNameCorrelation[[2]]
  print(sentrixID)
  sentrixPos <- str_remove_all(caseNameCorrelation[[3]], ".csv")
  print(sentrixPos)
  SUHMatch <- IDATSampleSheet %>% filter(Sentrix_ID == sentrixID, 
                                         Sentrix_Position == sentrixPos) %>% 
    dplyr::select(SUH) %>% dplyr::pull(SUH)
  
  cnvSes <- read.csv(outputDir) %>% 
    dplyr::select(c("chrom", "loc.start", "loc.end", "seg.mean")) %>% 
    dplyr::mutate(chrom = str_remove_all(chrom, "chr")) %>% 
    filter(!(chrom == "X") & !(chrom == "Y")) %>% 
    mutate(chrom = as.numeric(chrom)) %>% arrange(chrom) %>% 
    mutate(CNVStatus = case_when(
      seg.mean > 0.3 ~ "Amplification",
      seg.mean < -0.3 ~ "Deletion",
      TRUE ~ "Normal"), type = "SeSAMe") 
  cnvSes$chrom <- as.numeric(str_remove_all(cnvSes$chrom, pattern = "chr"))
  cnvSes <- cnvSes %>% filter(!(is.na(chrom))) %>% arrange(chrom)
  
  # SNP Processing
  chasFiltered <- chasFile[str_detect(chasFile$File, SUHMatch), ]
  chasFiltered <-  tidyr::separate_wider_delim(tidyr::separate_wider_delim(chasFiltered, `Full Location`, names = c("chromosome", "loc"), delim = ":"), `loc`,
                                               names = c("loc.start", "loc.end"), delim = "-") %>% dplyr::select(-c("chromosome")) %>% dplyr::rename(CNVStatus = "Type", seg.mean = "Mean Log2Ratio", chrom = "Chromosome") %>%
    mutate(
      type = "SNP", 
      CNVStatus = case_when(CNVStatus == "Loss" ~ "Deletion", 
                            CNVStatus == "Gain" ~ "Amplification", 
                            TRUE ~ "Normal")
    ) %>% dplyr::select(-c("File"))
  fasst2filtered <- fasst2file[str_detect(fasst2file$Sample, SUHMatch), ] %>%
    dplyr::select(c("Probe Median", "Event", "Chromosome Region"))
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
  
  rbind(chasFiltered, fasst2filtered, cnvSes) %>% arrange(chrom) %>% dplyr::filter(!(is.na(chrom)))
} # Sesame.
lmProcessingCon <- function(outputDir){
  IDATSampleSheet <- read.csv("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
  chasFile <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_ChAS_CNVs.xlsx") %>% dplyr::select(File, `Mean Log2Ratio`, Chromosome, Type, `Full Location`)
  fasst2file <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_FASST2_CNVs.xlsx")
  
  caseNameCorrelation <- str_split(str_split(outputDir, "/")[[1]][7], pattern = "_")[[1]]
  sentrixID <- caseNameCorrelation[[1]]
  print(sentrixID)
  sentrixPos <- str_remove_all(caseNameCorrelation[[2]], ".csv")
  print(sentrixPos)
  SUHMatch <- IDATSampleSheet %>% filter(Sentrix_ID == sentrixID, 
                                         Sentrix_Position == sentrixPos) %>% 
    dplyr::select(SUH) %>% dplyr::pull(SUH)
  
  cnvCon <- read.csv(outputDir) %>% 
    dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
  cnvCon$chrom <- as.numeric(str_remove_all(cnvCon$chrom, pattern = "chr"))
  cnvCon <- cnvCon %>% filter(!(is.na(chrom))) %>% arrange(chrom)
  
  # SNP Processing
  chasFiltered <- chasFile[str_detect(chasFile$File, SUHMatch), ]
  chasFiltered <-  tidyr::separate_wider_delim(tidyr::separate_wider_delim(chasFiltered, `Full Location`, names = c("chromosome", "loc"), delim = ":"), `loc`,
                                               names = c("loc.start", "loc.end"), delim = "-") %>% dplyr::select(-c("chromosome")) %>% dplyr::rename(CNVStatus = "Type", seg.mean = "Mean Log2Ratio", chrom = "Chromosome") %>%
    mutate(
      type = "SNP", 
      CNVStatus = case_when(CNVStatus == "Loss" ~ "Deletion", 
                            CNVStatus == "Gain" ~ "Amplification", 
                            TRUE ~ "Normal")
    ) %>% dplyr::select(-c("File"))
  fasst2filtered <- fasst2file[str_detect(fasst2file$Sample, SUHMatch), ] %>%
    dplyr::select(c("Probe Median", "Event", "Chromosome Region"))
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
  
  rbind(chasFiltered, fasst2filtered, cnvCon) %>% arrange(chrom) %>% dplyr::filter(!(is.na(chrom)))
} # Conumee.
labLMSProc <- function(STTq, Technology, binSize){
  basePath <- getwd()
  
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(paste0(basePath,"/LabData/LMS_SNP_EPIC_array_data/correlative.csv")) %>% filter(STT == STTq)
  
  cnvMatch <- read.csv(paste0(basePath,"/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV/", 
                              correlationSheet$CNV_Label[1], "_events.csv"))
  
  # methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
  #paste0(correlationSheet$Sentrix_ID[1],"_", 
  #     correlationSheet$Sentrix_Position[1], "/"), 
  # "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0(basePath,"/Outputs/Conumee/LabLMS/", binSize,"/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                                 correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0(basePath,"/Outputs/SeSAMe/LabLMS/", binSize,"/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
                                                                                               correlationSheet$Sentrix_Position[1],".csv")))
  
  cnvMatch <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(cnvMatch, `Chromosome.Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
                                          names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") %>% dplyr::rename(CNVStatus = "Event") %>%
    mutate(CNVStatus = case_when(
      CNVStatus == "CN Loss" ~ "Deletion", 
      CNVStatus == "CN Gain" ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% dplyr::rename(seg.mean = "Probe.Median") %>% 
    dplyr::select(c(chrom, loc.start, loc.end, CNVStatus, seg.mean, type)) %>% mutate(loc.start = as.numeric(gsub(",","",loc.start)), loc.end = as.numeric(gsub(",","",loc.end)))
  cnvMatch$chrom <- as.numeric(str_replace_all(cnvMatch$chrom, "chr", ""))
  cnvMatch <- cnvMatch %>% filter(!is.na(chrom))
  
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.2 ~ "Deletion", 
        seg.mean >= 0.2 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- rbind(methylMatch, cnvMatch) %>% mutate(Gene = as.character(0)) %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- rbind(conumeeMatch, cnvMatch) %>% mutate(Gene = as.character(0)) %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                 seg.mean < -0.3 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- rbind(sesameMatch, cnvMatch) %>% mutate(Gene = as.character(0)) %>% arrange(chrom)
  }
  combinedSet
}
#grangesSesame.R
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

#grangesMethyl.R
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
# Lab data.
labLMSProc <- function(STTq, Technology){
  basePath <- getwd()
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(paste0(basePath, "/LabData/LMS_SNP_EPIC_array_data/correlative.csv")) %>% filter(STT == STTq)
  
  cnvMatch <- read.csv(paste0(basePath,"/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV/", 
                              correlationSheet$CNV_Label[1], "_events.csv"))
  
  methylMatch <- read.csv(paste0(basePath,"/Outputs/MethylMaster/LabLMS/", 
                                 paste0(correlationSheet$Sentrix_ID[1],"_", 
                                        correlationSheet$Sentrix_Position[1], "/"), 
                                 "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0(basePath,"/Outputs/Conumee/LabLMS/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                    correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0(basePath,"/Outputs/SeSAMe/LabLMS/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
                                                                                  correlationSheet$Sentrix_Position[1],".csv")))
  
  cnvMatch <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(cnvMatch, `Chromosome.Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
                                          names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") %>% dplyr::rename(CNVStatus = "Event") %>%
    mutate(CNVStatus = case_when(
      CNVStatus == "CN Loss" ~ "Deletion", 
      CNVStatus == "CN Gain" ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% dplyr::rename(seg.mean = "Probe.Median") %>% 
    dplyr::select(c(chrom, loc.start, loc.end, CNVStatus, seg.mean, type)) %>% mutate(loc.start = as.numeric(gsub(",","",loc.start)), loc.end = as.numeric(gsub(",","",loc.end)))
  cnvMatch$chrom <- as.numeric(str_replace_all(cnvMatch$chrom, "chr", ""))
  cnvMatch <- cnvMatch %>% filter(!is.na(chrom))
  
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.2 ~ "Deletion", 
        seg.mean >= 0.2 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- rbind(methylMatch, cnvMatch) %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- rbind(conumeeMatch, cnvMatch) %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                 seg.mean < -0.3 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- rbind(sesameMatch, cnvMatch) %>% arrange(chrom)
  }
  
  
  
  combinedSet
}
labNmrlProc <- function(Sentrix, Technology, binSize){
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(
    "~/Work/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv") %>% filter(STT == STTq)
  
  # cnvMatch <- read.csv(paste0("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV/", 
  # correlationSheet$CNV_Label[1], "_events.csv"))
  
  # methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
  #paste0(correlationSheet$Sentrix_ID[1],"_", 
  #     correlationSheet$Sentrix_Position[1], "/"), 
  # "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0("~/Work/Analysis/Outputs/Conumee/LabNormals/", binSize,"/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                                     correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0("~/Work/Analysis/Outputs/SeSAMe/Normals/", binSize,"/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
                                                                                                correlationSheet$Sentrix_Position[1],".csv")))
  
  # cnvMatch <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(cnvMatch, `Chromosome.Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
  #                                        names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") %>% dplyr::rename(CNVStatus = "Event") %>%
  # mutate(CNVStatus = case_when(
  #   CNVStatus == "CN Loss" ~ "Deletion", 
  #   CNVStatus == "CN Gain" ~ "Amplification", 
  #   TRUE ~ "Normal"
  # )) %>% dplyr::rename(seg.mean = "Probe.Median") %>% 
  #  dplyr::select(c(chrom, loc.start, loc.end, CNVStatus, seg.mean, type)) %>% mutate(loc.start = as.numeric(gsub(",","",loc.start)), loc.end = as.numeric(gsub(",","",loc.end)))
  #cnvMatch$chrom <- as.numeric(str_replace_all(cnvMatch$chrom, "chr", ""))
  #cnvMatch <- cnvMatch %>% filter(!is.na(chrom))
  
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.2 ~ "Deletion", 
        seg.mean >= 0.2 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- methylMatch %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- conumeeMatch %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                 seg.mean < -0.3 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- sesameMatch %>% arrange(chrom)
  }
  
  
  
  combinedSet
}

#grangesConumee.R
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