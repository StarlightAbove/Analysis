### ============================================================
### 6. GENE CONCORDANCE
### ============================================================
# Per-oncogene log2-ratio agreement between each caller and the SNP array.
# Writes: quarterCutoff/Gene_concordance/<GENE>.csv, geneConcordance.csv

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and quarterCutoff/), since all data paths
# are built from getwd().
source("00_setup.R")
source("01_functions.R")
source("02_load_data.R")

# For a panel of oncogenes relevant to LMS/LM, compare each caller's
# log2 ratio at the gene to the SNP-array log2 ratio at the same gene,
# across bin sizes, via heatmaps and Pearson/Spearman correlations.

#### 6a. LMS ----
ap10000 <- NULL
ap50000 <- NULL
ap1e05 <- NULL
ap1e06 <- NULL

Gene <- c("MYC", "MYOCD", "CCNE1", "CDKN2A", "PTEN", "RB1", "TP53") 
bins <- c(10000, 50000, 1e+05, 1e+06)

for(stt in LMS_cases){
  ap10000 <- rbind(ap10000, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 10000), 
                                                                   labLMSProc(stt, "Conumee", 10000), 
                                                                   labLMSProc(stt, "Sesame", 10000)), 
                                           case = stt))
  ap50000 <- rbind(ap50000, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 50000), 
                                                                   labLMSProc(stt, "Conumee", 50000), 
                                                                   labLMSProc(stt, "Sesame", 50000)), 
                                           case = stt))
  ap1e05 <- rbind(ap1e05, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 1e+05), 
                                                                 labLMSProc(stt, "Conumee", 1e+05), 
                                                                 labLMSProc(stt, "Sesame", 1e+05)), 
                                         case = stt))
  ap1e06 <- rbind(ap1e06, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 1e+06), 
                                                                 labLMSProc(stt, "Conumee", 1e+06), 
                                                                 labLMSProc(stt, "Sesame", 1e+06)), 
                                         case = stt))
}
techs <- c("Gene_MMasteR", "Gene_SNP", "Gene_Conumee", "Gene_SeSAMe")

ap10000 <- ap10000 %>% dplyr::rename(log2ratio = seg.mean)
ap50000 <- ap50000 %>% dplyr::rename(log2ratio = seg.mean)
ap1e05 <- ap1e05 %>% dplyr::rename(log2ratio = seg.mean)
ap1e06 <- ap1e06 %>% dplyr::rename(log2ratio = seg.mean)

ggplot(ap10000, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 10000") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(ap50000, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 50000") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(ap1e05, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 1e+05") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(ap1e06, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 1e+06") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()


#### 6b. Save raw per-gene data ----
ap10000 <- ap10000 %>% dplyr::mutate(bin = 10000)
ap50000 <- ap50000 %>% dplyr::mutate(bin = 50000)
ap1e05 <- ap1e05 %>% dplyr::mutate(bin = 1e+05)
ap1e06 <- ap1e06 %>% dplyr::mutate(bin = 1e+06)
ap <- rbind(ap10000, ap50000, ap1e05, ap1e06) %>% dplyr::select(-c("chrom")) %>%
  dplyr::mutate(width = loc.end - loc.start)
aps <- split(ap, ap$Gene)

for (i in seq_along(aps)) {
  file_name <- paste0(names(aps)[i], ".csv")
  write.csv(aps[[i]], file = paste0(getwd(),"/quarterCutoff/Gene_concordance/", file_name), row.names = FALSE)
}
aps[[Gene[1]]] <- split(aps[[Gene[1]]], aps[[Gene[1]]]$bin)
aps[[Gene[2]]] <- split(aps[[Gene[2]]], aps[[Gene[2]]]$bin)
aps[[Gene[3]]] <- split(aps[[Gene[3]]], aps[[Gene[3]]]$bin)
aps[[Gene[4]]] <- split(aps[[Gene[4]]], aps[[Gene[4]]]$bin)
aps[[Gene[5]]] <- split(aps[[Gene[5]]], aps[[Gene[5]]]$bin)
aps[[Gene[6]]] <- split(aps[[Gene[6]]], aps[[Gene[6]]]$bin)
aps[[Gene[7]]] <- split(aps[[Gene[7]]], aps[[Gene[7]]]$bin)

#### 6c. Concordance correlations (caller vs. SNP truth, per gene/bin) ----
# NOTE: this reassigns the global `bins` to character values; it stays
# reassigned for the rest of the script (see the Accuracy section below).
bins <- c("50000", "10000", "1e+05", "1e+06")

snp_data <- ap %>%
  filter(type == "Gene_SNP") %>%
  select(bin, case, Gene, log2ratio) %>%
  rename(log2ratio_SNP = log2ratio)

# Step 2: Non-SNP data — average out duplicates per bin + case + Gene + type
non_snp_data <- ap %>%
  filter(type != "Gene_SNP") %>%
  group_by(bin, case, Gene, type) %>%
  summarise(
    n_averaged = n(),
    log2ratio  = mean(log2ratio, na.rm = TRUE),
    .groups    = "drop"
  )

# Step 3: Join to SNP by bin + case + Gene
paired <- non_snp_data %>%
  inner_join(snp_data, by = c("bin", "case", "Gene"))

# Step 4: Correlate each type vs SNP, faceted by bin + Gene
correlations <- paired %>%
  group_by(bin, Gene, type) %>%
  summarise(
    n          = n(),
    pearson_r  = cor(log2ratio, log2ratio_SNP, method = "pearson"),
    spearman_r = cor(log2ratio, log2ratio_SNP, method = "spearman"),
    .groups    = "drop"
  )

# Step 5: Human readable formatting.
output.df.pearson <- correlations %>%
  mutate(Tech = case_when(
    type == "Gene_MMasteR"  ~ "MethylMasteR",
    type == "Gene_Conumee"  ~ "Conumee",
    type == "Gene_SeSAMe"   ~ "Sesame",
    TRUE ~ type  # fallback: keep as-is
  )) %>%
  select(Tech, Bin = bin, Gene, pearson_r) %>%
  pivot_wider(
    names_from  = Gene,
    values_from = pearson_r
  ) %>%
  arrange(Tech, Bin)

output.df.spearman <- correlations %>%
  mutate(Tech = case_when(
    type == "Gene_MMasteR"  ~ "MethylMasteR",
    type == "Gene_Conumee"  ~ "Conumee",
    type == "Gene_SeSAMe"   ~ "Sesame",
    TRUE ~ type  # fallback: keep as-is
  )) %>%
  select(Tech, Bin = bin, Gene, spearman_r) %>%
  pivot_wider(
    names_from  = Gene,
    values_from = spearman_r
  ) %>%
  arrange(Tech, Bin)

make_wide <- function(metric_col, metric_name) {
  correlations %>%
    mutate(Tech = case_when(
      type == "Gene_MMasteR" ~ "MethylMasteR",
      type == "Gene_Conumee" ~ "Conumee",
      type == "Gene_SeSAMe"  ~ "Sesame",
      TRUE ~ type
    ),
    Bin = factor(bin, levels = c(10000, 50000, 100000, 1000000))
    ) %>%
    select(Tech, Bin, Gene, value = {{ metric_col }}) %>%
    mutate(Metric = metric_name) %>%
    pivot_wider(
      names_from  = Gene,
      values_from = value
    ) %>%
    arrange(Tech, Bin)
}

result_pearson  <- make_wide(pearson_r,  "Pearson")
result_spearman <- make_wide(spearman_r, "Spearman")

result_long_alt <- bind_rows(result_pearson, result_spearman) %>%
  arrange(Tech, Bin, Metric) %>%
  mutate(
    Tech_Bin = paste(Tech, Bin, sep = " | "),
    Tech_Bin = factor(Tech_Bin, levels = unique(Tech_Bin))
  ) %>%
  select(-Tech, -Bin) %>%
  pivot_longer(
    cols      = -c(Tech_Bin, Metric),
    names_to  = "Gene",
    values_to = "Correlation"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("Pearson", "Spearman")))

# Plot: Gene on Y, Tech_Bin faceted as columns, Metric as x sub-groups
ggplot(result_long_alt, aes(x = Metric, y = Gene, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Correlation, 2)), size = 2.5) +
  scale_fill_gradient2(
    low      = "red",
    mid      = "pink",
    high     = "green",
    limits   = c(min(result_long_alt$Correlation, na.rm = TRUE), 
                 max(result_long_alt$Correlation, na.rm = TRUE)),
    midpoint = mean(result_long_alt$Correlation, na.rm = TRUE),
    name     = "Correlation"
  ) +
  facet_grid(
    cols   = vars(Tech_Bin),
    switch = "x"                  # puts Tech|Bin labels at the bottom
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y       = element_text(size = 9),
    strip.text.x      = element_text(angle = 90, hjust = 0, size = 8),  # rotate facet labels
    panel.spacing     = unit(0.1, "lines"),  # tighten columns
    panel.grid        = element_blank()
  ) +
  labs(
    title = "Correlation Heatmap: Non-SNP Methods vs SNP",
    x     = NULL,
    y     = "Gene"
  )

output.df.pearson <- output.df.pearson %>% dplyr::mutate(Metric = "Pearson")
output.df.spearman <- output.df.spearman %>% dplyr::mutate(Metric = "Spearman")
output.df <- rbind(output.df.pearson, output.df.spearman)
write.csv(output.df, "quarterCutoff/Gene_concordance/geneConcordance.csv")
