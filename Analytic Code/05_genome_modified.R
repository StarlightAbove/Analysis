### ============================================================
### 5. GENOME MODIFIED
### ============================================================
# % of genome altered per case / bin / caller for LMS and LM, vs. SNP baseline.
# Writes: results/Genome_modified/genome_modifiedLMS.csv, genome_modifiedLM.csv

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and results/), since all data paths
# are built from getwd().
source(paste0(getwd(), "/Analytic Code/00_setup.R"))
source(paste0(getwd(), "/Analytic Code/01_functions.R"))
source(paste0(getwd(), "/Analytic Code/02_load_data.R"))

# % of genome altered (see GenomeModified()) per case/bin/caller, plus the
# corresponding SNP-array baseline pulled from the ChAS design spreadsheet.

#### 5a. LMS ----
techs <- c("MethylMaster", "Sesame", "Conumee")
bins <- c(50000, 10000, 1e+05, 1e+06)
geneMod <- NULL
for(cases in LMS_cases){
  for(bin in bins){
    for(tech in techs){
      result <- data.frame(
        category = paste0(cases, "-", bin, "-", tech),
        val = GenomeModified(labLMSProc(STTq = cases, Technology = tech, binSize = bin) )
      )
      
      geneMod <- rbind(geneMod, result)
    }
  }
}

chasFile <- read_excel("LabData/LMS_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/design_13LMS_CNVs_other_info_01Feb2026.xlsx") %>%
  dplyr::select(c(STT, "% Genome Changed")) %>%
  dplyr::rename(Case = STT, val = "% Genome Changed") %>%
  dplyr::mutate(Bin = "DEF", Tech = "SNP")


geneMod <- geneMod %>%
  separate(
    col = category,
    into = c("Case", "Bin", "Tech"),
    sep = "-"
  )

geneMod <- rbind(geneMod, chasFile)
geneMod <- geneMod %>% arrange(Case)
write.csv(geneMod, file = paste0(getwd(), "/results/Genome_modified/genome_modifiedLMS.csv"))
df <- read.csv("results/Genome_modified/genome_modifiedLMS.csv") %>% dplyr::select(-c("X"))
snp_data <- df %>% filter(Tech == "SNP" & Bin == "DEF")
main_data <- df %>% filter(Bin != "DEF")
main_data$Bin <- factor(main_data$Bin, levels = c("10000", "50000", "1e+05", "1e+06", "1e+07"))
plt <- ggplot() +
  # SNP baseline: horizontal reference line per Case
  geom_hline(
    data = snp_data,
    aes(yintercept = val, linetype = "SNP baseline"),
    colour = "grey40",
    linewidth = 0.7
  ) +
  # Main lines for the three tools
  geom_line(
    data = main_data,
    aes(x = Bin, y = val, colour = Tech, group = Tech),
    linewidth = 0.8
  ) +
  geom_point(
    data = main_data,
    aes(x = Bin, y = val, colour = Tech),
    size = 2
  ) +
  # Facet by Case
  facet_wrap(~ Case, ncol = 4) +
  # Scales and labels
  scale_colour_brewer(palette = "Set1", name = "Tool") +
  scale_linetype_manual(
    name   = NULL,
    values = c("SNP baseline" = "dashed")
  ) +
  labs(
    title = "Genome Modified Across Bins and Cases",
    subtitle = "Dashed line = per-case SNP baseline value",
    x     = "Bin size",
    y     = "Value"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

#### 5b. LM ----
techs <- c("MethylMaster", "Sesame", "Conumee")
bins <- c(50000, 10000, 1e+05, 1e+06)
geneModLM <- NULL
for(tec in techs){
  for(b in bins){
    for(cases in LM_cases){
      result <- data.frame(
        category = paste0(cases, "-", b, "-", tec),
        val = GenomeModified(LMStt(STT = cases, bin = b, tech = tec) )
      )
      geneModLM <- rbind(geneModLM, result)
    }
  }
}

GenomeModifiedSNP <- NULL
for(cases in LM_cases){
  
  SNPDf <- read_delim(paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_LM_Probe_and_segment_level_data_01Feb2026/STT",
                      cases, "_Segment_level_data_01Feb2026.segment.txt")) %>% 
    dplyr::select(c("Chromosome", "StartPosition", "StopPosition" ,"Median Log2 Ratio")) %>%
    drop_na() %>%
    dplyr::rename(seg.mean = "Median Log2 Ratio") %>%
    dplyr::filter(seg.mean < -0.25 | seg.mean > 0.25) %>%
    dplyr::mutate(width = StopPosition - StartPosition) %>%
    dplyr::filter(width > 1e+06) #1Mb cut off
  
  widthSum <- sum(SNPDf$width)
  
  result <- data.frame(
    category = paste0(cases, "-DEF-SNP"),
    val = widthSum/hg19_total
  )
  GenomeModifiedSNP <- rbind(GenomeModifiedSNP, result)
}

snp_data <- GenomeModifiedSNP %>%
  separate(
    col = category,
    into = c("Case", "Bin", "Tech"),
    sep = "-"
  )

main_data <- geneModLM %>%
  separate(
    col = category,
    into = c("Case", "Bin", "Tech"),
    sep = "-"
  )
main_data$Bin <- factor(main_data$Bin, levels = c("10000", "50000", "1e+05", "1e+06"))

ggplot() +
  # SNP baseline: horizontal reference line per Case
  geom_hline(
    data = snp_data,
    aes(yintercept = val, linetype = "SNP baseline"),
    colour = "grey40",
    linewidth = 0.7
  ) +
  # Main lines for the three tools
  geom_line(
    data = main_data,
    aes(x = Bin, y = val, colour = Tech, group = Tech),
    linewidth = 0.8
  ) +
  geom_point(
    data = main_data,
    aes(x = Bin, y = val, colour = Tech),
    size = 2
  ) +
  # Facet by Case
  facet_wrap(~ Case, ncol = 4) +
  # Scales and labels
  scale_colour_brewer(palette = "Set1", name = "Tool") +
  scale_linetype_manual(
    name   = NULL,
    values = c("SNP baseline" = "dashed")
  ) +
  labs(
    title = "Genome Modified Across Bins and Cases",
    subtitle = "Dashed line = per-case SNP baseline value",
    x     = "Bin size",
    y     = "Value"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

LMdf <- rbind(snp_data, main_data)
write.csv(LMdf, file = paste0(getwd(), "/results/Genome_modified/genome_modifiedLM.csv"))
