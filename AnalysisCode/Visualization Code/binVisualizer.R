labLMSProc <- function(STTq, Technology, binSize){
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv") %>% filter(STT == STTq)
  
  cnvMatch <- read.csv(paste0("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV/", 
                              correlationSheet$CNV_Label[1], "_events.csv"))
  
  # methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
  #paste0(correlationSheet$Sentrix_ID[1],"_", 
  #     correlationSheet$Sentrix_Position[1], "/"), 
  # "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0("~/Work/Analysis/Outputs/Conumee/LabLMS/", binSize,"/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                                 correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0("~/Work/Analysis/Outputs/SeSAMe/LabLMS/", binSize,"/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
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
binLength <-  c(10000,25000,50000,75000,100000, 5e+05, 1e+06)


lablmsCodes <- read.csv("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
LabNmrlsCodes <- read.csv(
  "~/Work/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")$Basename

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
  
  # Vertical chromosome boundaries
  chr_boundaries <- chr_lengths %>%
    mutate(x = chr_start) %>%
    select(chrom, x)
  
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
  ggsave(plot = plt1, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/25000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt2, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/50000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt3, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/75000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt4, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/100000",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt5, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/5e+05",
         width = 1920, height = 1080, units = "px")
  ggsave(plot = plt6, filename = paste0(toString(nmls), ".png"), 
         path = "~/Work/Analysis/Statistics/Normals/byBin/1e+06",
         width = 1920, height = 1080, units = "px")
  
}

# Plotting for TCGA LMS.

# Plotting for LM.
