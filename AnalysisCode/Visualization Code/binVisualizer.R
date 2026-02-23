library(patchwork)
library(cowplot)
library(gridGraphics)
library(gridExtra)
library(ggplotify)
labLMSProc <- function(STTq, Technology, binSize){
  
  
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv") %>% filter(STT == STTq)
  
  cnvMatch <- read_delim(paste0("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_LMS_Probe_and_segment_level_data_01Feb2026/STT", 
                              STTq, "_Recentered_Segment_level_data_01Feb2026.segment.txt"))
  
  # methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
  #paste0(correlationSheet$Sentrix_ID[1],"_", 
  #     correlationSheet$Sentrix_Position[1], "/"), 
  # "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0("~/Work/Analysis/Outputs/Conumee/LabLMS/", binSize,"/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                                 correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0("~/Work/Analysis/Outputs/SeSAMe/LabLMS/", binSize,"/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
                                                                                               correlationSheet$Sentrix_Position[1],".csv")))
  
  methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", binSize,"/", 
                                 correlationSheet$Sentrix_ID[1],"_", correlationSheet$Sentrix_Position[1],"/autocorrected_regions.csv"))
  
  cnvMatch <- cnvMatch %>% dplyr::filter(!(Type == "LOH")) %>% dplyr::filter(Chromosome != 24 & Chromosome != 25) %>% dplyr::select("Chromosome", "StartPosition", "StopPosition", "Median Log2 Ratio") %>% 
    dplyr::rename(chrom = "Chromosome", loc.start = "StartPosition", loc.end = "StopPosition", seg.mean = "Median Log2 Ratio") %>%
    dplyr::mutate(seg.mean = as.numeric(seg.mean)) %>%
    dplyr::mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    ), type = "SNP")
  
  print(cnvMatch)
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
    combinedSet <- rbind(sesameMatch, cnvMatch) %>% mutate(Gene = as.character(0), seg.mean = as.numeric(seg.mean)) %>% arrange(chrom)
  }
  combinedSet
}
Gene <- c("MYC", "MYOCD", "CCNE1", "CDKN2A", "PTEN", "RB1", "TP53") 
geneAnno <- function(Gene, db = NULL){
  reference <- geneGen(Gene = Gene, db = db)
  print(reference[[1]])
  pheatmap_ggplot <- as.ggplot(reference[[1]]$gtable)
  reference2 <- reference[[2]]
  db <- unique(rbind(db, reference2))
  img2 <- plot_cnv_segments(db)
  return(list((img2 + pheatmap_ggplot), db))
}
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


library(tidyverse)
library(reshape2)

binLength <- c(10000,25000,50000,75000,1e+05, 5e+05, 1e+06)


lablmsCodes <- read.csv("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
LabNmrlsCodes <- read.csv(
  "~/Work/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")$Basename
plot_cnv_segments(rbind(labLMSProc(lablmsCodes[1], "Conumee", 10000), 
                        labLMSProc(lablmsCodes[1], "Sesame", 10000)))
rbind(labLMSProc(lablmsCodes[1], "Conumee", 10000), 
      labLMSProc(lablmsCodes[1], "Sesame", 10000))
labLMSPlots <- c()
labNmrlPlots <- c()

for(llC in lablmsCodes){
  plt <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 10000), 
                                 labLMSProc(llC, "Sesame", 10000)))
  plt1 <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 25000), 
                                 labLMSProc(llC, "Sesame", 25000)))
  plt2 <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 50000), 
                                 labLMSProc(llC, "Sesame", 50000)))
  plt3 <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 75000), 
                                 labLMSProc(llC, "Sesame", 75000)))
  plt4 <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 100000), 
                                 labLMSProc(llC, "Sesame", 100000)))
  plt5 <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 5e+05), 
                                 labLMSProc(llC, "Sesame", 5e+05)))
  plt6 <- plot_cnv_segments(rbind(labLMSProc(llC, "Conumee", 1e+06), 
                                 labLMSProc(llC, "Sesame", 1e+06)))
  
  ggsave(plot = plt, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/10000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt1, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/25000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt2, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/50000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt3, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/75000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt4, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/100000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt5, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/5e+05",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt6, filename = paste0(toString(llC), ".png"), 
         path = "~/Work/Analysis/Statistics/LabLMS/byBin/1e+06",
         width = 1920, height = 1080, units = "px")
}

plot_cnvns(rbind(labNmrlProc(LabNmrlsCodes[1],"Sesame","10000"), 
                 labNmrlProc(LabNmrlsCodes[1],"Conumee","10000")))

# Plotting for lab normals.
plot_cnvns <- function(df){
  df <- df %>%
    mutate(
      type_group = ifelse(type == "SNP", "SNP", "non-SNP")
    ) %>%
    arrange(chrom, loc.start)
  
  # Calculate chromosome cumulative positions
  chr_lengths <- df %>%
    group_by(chrom) %>%
    summarize(chr_len = max(loc.end), .groups = "drop") %>%
    arrange(chrom) %>%
    mutate(chr_start = lag(cumsum(as.numeric(chr_len)), default = 0)) %>%
    mutate(chr_mid = chr_start + chr_len / 2)
  
  # Join to get cumulative start and end positions
  df <- df %>%
    left_join(chr_lengths, by = "chrom") %>%
    mutate(
      start_cum = loc.start + chr_start,
      end_cum = loc.end + chr_start
    ) %>%mutate(type = ifelse(type == "SNP", "SNP Array", type)) %>% mutate(type = as.factor(type))
  print(df)
  # Vertical chromosome boundaries
  chr_boundaries <- chr_lengths %>%
    mutate(x = chr_start) %>%
    dplyr::select(chrom, x)
  
  x_breaks <- chr_lengths$chr_mid
  x_labels <- paste0("chr", chr_lengths$chrom)
  print(x_labels)
  print(df)
  
  p <- ggplot(df, aes(x = start_cum, xend = end_cum, y = seg.mean, yend = seg.mean)) +
    geom_segment(aes(color = type), size = 0.7, alpha = 0.8) +
    geom_vline(data = chr_boundaries, aes(xintercept = x), color = "grey70", linetype = "dashed") +
    # scale_color_manual(values = c("Amplification" = "red", "Deletion" = "blue", "Normal" = "black")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("SeSAMe" = "red", "MethylMaster" = "blue", "Conumee" = "green", "SNP Array" = "black")) +
    labs(
      x = "Genomic Position (across chromosomes)",
      y = "Segment Mean (log2 ratio)",
      title = "CNV Segments Across Genome"
    ) + geom_hline(yintercept = -0.2, linetype = "dotted", color = "black") + 
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "black") + 
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )
  
  return(p)
}
for(nmls in LabNmrlsCodes){
  plt <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","10000"), 
                   labNmrlProc(nmls,"Conumee","10000")))
  plt2 <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","25000"), 
                          labNmrlProc(nmls,"Conumee","25000")))
  plt3 <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","50000"), 
                          labNmrlProc(nmls,"Conumee","50000")))
  plt4 <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","75000"), 
                          labNmrlProc(nmls,"Conumee","75000")))
  plt5 <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","1e+05"), 
                          labNmrlProc(nmls,"Conumee","1e+05")))
  plt6 <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","5e+05"), 
                          labNmrlProc(nmls,"Conumee","5e+05")))
  plt7 <- plot_cnvns(rbind(labNmrlProc(nmls,"Sesame","1e+06"), 
                          labNmrlProc(nmls,"Conumee","1e+06")))
  ggsave(plot = plt, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/10000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt2, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/25000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt3, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/50000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt4, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/75000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt5, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/100000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt6, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/5e+05",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt7, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/1e+06",
         width = 1920, height = 1080, units = "px")
  
}

## Plotting for LM.
binLength <- c(10000,25000,50000,75000,1e+05, 5e+05, 1e+06)

# Getting paired paths.
corr <- read.csv("~/Work/Analysis/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv") %>%
        dplyr::select(c(Sentrix_ID, Sentrix_Position))
corr <-paste0(corr$Sentrix_ID, "_", corr$Sentrix_Position, ".csv")
basePthSes <- "./Outputs/SeSAMe/LM/bins/"
basePthCon <- "./Outputs/Conumee/LMData/bins/"
basePthSes <- paste0(basePthSes, binLength, "/segments_")
basePthCon <- paste0(basePthCon, binLength, "/")
pathsSes <- c()
pathsCon <- c()
for(bps in basePthSes){
  pathsSes <- c(pathsSes, paste0(bps, corr))
}
for(bps in basePthCon){
  pathsCon <- c(pathsCon, paste0(bps, corr))
}
paired <- Map(list, pathsSes, pathsCon)
names(paired) <- c(1:63)
graphs <- c()
for(pr in paired){
  graphs <- c(graphs, 
              plot_cnv_segments(rbind(lmProcessingCon(pr[[2]]), 
                                      lmProcessingSes(pr[[1]]))))
}
for(i in 1:9){
  ggsave(plot = graphs[[1*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/10000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = graphs[[2*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/25000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = graphs[[3*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/50000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = graphs[[4*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/75000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = graphs[[5*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/100000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = graphs[[6*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/5e+05",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = graphs[[7*i]], filename = paste0(corr[i], ".png"), 
         path = "~/Work/Analysis/Statistics/LM/byBin/1e+06",
         width = 1920, height = 1080, units = "px")
}

# Bin differentials.
case1 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "10000"))
case2 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "25000"))
case3 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "50000"))
case4 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "75000"))
case5 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "1e+05"))
case6 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "5e+05"))
case7 <- accuracyModel(labLMSProc(lablmsCodes[1], "Conumee", "1e+06"))
accuracies <- data.frame(chrom = as.numeric(case1[[1]]$Chromosome), 
                         "10000" = case1[[1]]$Accuracy,
                         "25000" = case2[[1]]$Accuracy,
                         "50000" = case3[[1]]$Accuracy,
                         "75000" = case4[[1]]$Accuracy,
                         "1e+05" = case5[[1]]$Accuracy,
                         "5e+05" = case6[[1]]$Accuracy,
                         "1e+06" = case7[[1]]$Accuracy)%>% 
  dplyr::arrange(chrom) %>% 
  pivot_longer(
    cols = starts_with("X"), 
    names_to = "BinSize", 
    values_to = "Accuracies"  
  )
ggplot(accuracies, aes(group = chrom, y = Accuracies, x = chrom)) +
  geom_boxplot() + facet_wrap(~BinSize, scale="free") + geom_point() +
  labs(title = "STT 9202 Accuracy Across Different Bins")


case1 <- fpCheck(lmProcessingCon(pathsCon[1]))
case2 <- fpCheck(lmProcessingCon(pathsCon[9]))
case3 <- fpCheck(lmProcessingCon(pathsCon[18]))
case4 <- fpCheck(lmProcessingCon(pathsCon[27]))
case5 <- fpCheck(lmProcessingCon(pathsCon[36]))
case6 <- fpCheck(lmProcessingCon(pathsCon[45]))
case7 <- fpCheck(lmProcessingCon(pathsCon[54]))
accuracies <- data.frame(chrom = as.numeric(case1[[1]]$Chromosome), 
                         "10000" = case1[[1]]$Accuracy,
                         "25000" = case2[[1]]$Accuracy,
                         "50000" = case3[[1]]$Accuracy,
                         "75000" = case4[[1]]$Accuracy,
                         "1e+05" = case5[[1]]$Accuracy,
                         "5e+05" = case6[[1]]$Accuracy,
                         "1e+06" = case7[[1]]$Accuracy)%>% 
  dplyr::arrange(chrom) %>% 
  pivot_longer(
    cols = starts_with("X"), 
    names_to = "BinSize", 
    values_to = "Accuracies"  
  )
ggplot(accuracies, aes(group = chrom, y = Accuracies, x = chrom)) +
  geom_boxplot() + facet_wrap(~BinSize, scale="free") + geom_point() +
  labs(title = "203219650162_R08C01 Accuracy Across Different Bins")
