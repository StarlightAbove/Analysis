### ============================================================
### 4. GENOMIC INDEX
### ============================================================
# Genomic Index per case / bin / caller for LMS and LM, plus bar plots.
# Writes: results/Genomic_Index/genomicIndexLMS.csv, genomicIndexLM.csv

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and results/), since all data paths
# are built from getwd().
source(paste0(getwd(), "/Analytic Code/00_setup.R"))
source(paste0(getwd(), "/Analytic Code/01_functions.R"))
source(paste0(getwd(), "/Analytic Code/02_load_data.R"))

# Genomic Index (see GenomicIndexIntermediateMatrix()/GenomicIndex()) computed
# per case, per bin size, per caller (incl. SNP truth at 1Mb), for both cohorts.

#### 4a. LMS ----
outputDf <- NULL
for(i in LMS_cases){
  print(i)
  for(j in bins){
    print(j)
    mm <- GenomicIndexIntermediateMatrix(labLMSProc(i, 
                                                    Technology = "MethylMaster", 
                                                    binSize = j), 
                                         stt = i, bin = j, LMSorLM = "LMS", 
                                         tech = "MethylMaster")
    cn <- GenomicIndexIntermediateMatrix(labLMSProc(i, 
                                                    Technology = "Sesame", 
                                                    binSize = j), 
                                         stt = i, bin = j, LMSorLM = "LMS", 
                                         tech = "Sesame")
    ss <- GenomicIndexIntermediateMatrix(labLMSProc(i, 
                                                    Technology = "Conumee", 
                                                    binSize = j), 
                                         stt = i, bin = j, LMSorLM = "LMS", 
                                         tech = "Conumee")
    outputDf <- unique(rbind(outputDf, mm, cn, ss))
  }
}

summ <- outputDf %>% group_by(Case, Bin, type) %>% summarize(
  count = n(),
  c_sq = (n())^2,
  chr_count = n_distinct(chrom)
) %>% dplyr::mutate(gi = c_sq/chr_count) %>% filter(!(Bin != 1e+06 & type == "SNP"))

write.csv(summ, paste0(getwd(), "/results/Genomic_Index/genomicIndexLMS.csv"))

#### 4b. LM ----
outputDf <- NULL
for(i in LM_cases){
  for(j in bins){
    mm <- GenomicIndexIntermediateMatrix(LMStt(STT = i, bin = j, tech = "MethylMaster"), 
                                         stt = i, bin = j, LMSorLM = "LM", 
                                         tech = "MethylMaster")
    cn <- GenomicIndexIntermediateMatrix(LMStt(STT = i, bin = j, tech = "Sesame"), 
                                         stt = i, bin = j, LMSorLM = "LM", 
                                         tech = "Sesame")
    ss <- GenomicIndexIntermediateMatrix(LMStt(STT = i, bin = j, tech = "Conumee"), 
                                         stt = i, bin = j, LMSorLM = "LM", 
                                         tech = "Conumee")
    outputDf <- unique(rbind(outputDf, mm, cn, ss))
  }
}

LMGI <- outputDf %>% group_by(Case, Bin, type) %>% summarize(
  count = n(),
  c_sq = (n())^2,
  chr_count = n_distinct(chrom)
) %>% dplyr::mutate(gi = c_sq/chr_count)

write.csv(LMGI, paste0(getwd(), "/results/Genomic_Index/genomicIndexLM.csv"))

#### 4c. Plotting ----
# Bar chart of 1 Mb Genomic Index by case, colored by caller. Case/tool
# combinations with no computed value (no row in the source CSV) are filled
# in as gi = 0 and drawn as a flat line at zero instead of a bar, so a
# genuine zero result is visually distinguishable from a missing one.
giplot <- function(df){
  p <- ggplot(df, aes(x = Case, y = gi, fill = type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = 0, ymax = 0, alpha = filled, color = type),
                  position = position_dodge(width = 0.8), width = 0.7,
                  linewidth = 0.6, show.legend = FALSE) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none") +
    labs(
      title = "Genomic Index by Case and Tool (1 Mb bins)",
      x     = "Case",
      y     = "Genomic Index",
      fill  = "Tool"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey92"),
      legend.position  = "bottom"
    )
}

df <- read.csv(paste0(getwd(), "/results/Genomic_Index/genomicIndexLM.csv")) %>%
  select(Case, Bin, type, gi) %>%
  filter(Bin == 1e+06) %>%
  mutate(Case = factor(Case)) %>%
  tidyr::complete(Case, type, fill = list(gi = 0)) %>%
  mutate(filled = gi == 0)
LMPlot <- giplot(df)

df <- read.csv(paste0(getwd(), "/results/Genomic_Index/genomicIndexLMS.csv")) %>%
  select(Case, Bin, type, gi) %>%
  filter(Bin == 1e+06) %>%
  mutate(Case = factor(Case)) %>%
  tidyr::complete(Case, type, fill = list(gi = 0)) %>%
  mutate(filled = gi == 0)
LMSPlot <- giplot(df)
