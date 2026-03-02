# Bioconductor first
library(GenomeInfoDb)
library(AnnotationDbi)
library(AnnotationHub)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rentrez)
library(GenomicRanges)
library(GenomicFeatures)
library(pheatmap)

# Then tidyverse and dplyr last so they win all masking conflicts
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(gridGraphics)
library(gridExtra)
library(ggplotify)
library(readxl)
library(RColorBrewer)
library(reshape2)

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
LMStt <- function(STT, bin, tech){
  
  corrSheet <- read.csv("~/Work/Analysis/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv") %>%
    dplyr::filter(STT == STT)
  sentrixID <- corrSheet$Sentrix_ID[1]
  sentrixPos <- corrSheet$Sentrix_Position[1]
  
  SNP <- read_delim(paste0("~/Work/Analysis/LabData/LM_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_LM_Probe_and_segment_level_data_01Feb2026/STT",
                           STT,"_Segment_level_data_01Feb2026.segment.txt"))
  SNP <- SNP %>% dplyr::filter(!(Type == "LOH")) %>% dplyr::filter(Chromosome != 24 & Chromosome != 25) %>% dplyr::select("Chromosome", "StartPosition", "StopPosition", "Median Log2 Ratio") %>% 
    dplyr::rename(chrom = "Chromosome", loc.start = "StartPosition", loc.end = "StopPosition", seg.mean = "Median Log2 Ratio") %>%
    dplyr::mutate(seg.mean = as.numeric(seg.mean)) %>%
    dplyr::mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    ), type = "SNP")
  
  
  # MethylMasteR
  cnvMethyl <- NULL;
  if(tech == "MethylMaster"){
    outputDirMethyl <- paste0("~/Work/Analysis/Outputs/MethylMaster/LM/",bin,"/",sentrixID,"_",sentrixPos,"/autocorrected_regions.csv")
    cnvMethyl <- read.csv(outputDirMethyl) %>%
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
  }
  
  
  # SeSAMe
  sesameOutput <- NULL
  if(tech == "Sesame"){
    outputDirSesame <- paste0("~/Work/Analysis/Outputs/SeSAMe/LM/bins/",bin,"/","segments_",sentrixID, "_", sentrixPos, ".csv")
    sesameOutput <- read.csv(outputDirSesame) %>% dplyr::select(c("chrom", "loc.start", 
                                                                  "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                                    chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                                      !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                                        chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                                          CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                                seg.mean < -0.3 ~ "Deletion",
                                                                                                TRUE ~ "Normal"), type = "SeSAMe") 
  }
  
  
  # Conumee
  case <- NULL
  if(tech == "Conumee"){
    outputDirConumee <- paste0("~/Work/Analysis/Outputs/Conumee/LMData/bins/",bin,"/",sentrixID,"_",sentrixPos,".csv")
    case <- read.csv(outputDirConumee)
    case <- case %>% dplyr::select(-c("ID", "bstat")) %>% 
      mutate(CNVStatus = case_when(case$seg.mean > 0.2 ~ "Amplification",
                                   case$seg.mean < -0.2 ~ "Deletion",
                                   TRUE ~ "Normal"), 
             chrom = as.numeric(str_remove_all(chrom, "chr"))) %>% 
      dplyr::arrange(chrom) %>% mutate(type = "Conumee") %>%
      dplyr::select(-c("num.mark", "pval", "seg.median", "X"))
  }
  
  output <- rbind(SNP,case,sesameOutput,cnvMethyl)
  output$chrom <- as.character(output$chrom)
  output
}
labNmrlProc <- function(Sentrix, Technology, binSize){
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(
    "~/Work/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv") %>% filter(Basename == Sentrix)
  
  conumeeMatch <- read.csv(paste0("~/Work/Analysis/Outputs/Conumee/LabNormals/", binSize,"/", paste0(Sentrix,".csv")))
  
  sesameMatch <- read.csv(paste0("~/Work/Analysis/Outputs/SeSAMe/Normals/", binSize,"/", paste0("segments_",Sentrix,".csv")))
  
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

