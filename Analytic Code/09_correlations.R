### ============================================================
### 9. CORRELATIONS
### ============================================================
# Accuracy vs. Genomic Index / Genome Modified (bootstrap Spearman),
# and SNP vs. methylation agreement on both metrics.
#
# Requires outputs from: 04_genomic_index.R, 05_genome_modified.R, 07_accuracy.R
# Writes: results/Correlations/*.csv

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and results/), since all data paths
# are built from getwd().
source(paste0(getwd(), "/Analytic Code/00_setup.R"))
source(paste0(getwd(), "/Analytic Code/01_functions.R"))
source(paste0(getwd(), "/Analytic Code/02_load_data.R"))

# Bootstrap Spearman correlations between Accuracy and Genomic
# Index / Genome Modified, plus direct SNP-vs-methylation agreement
# on Genomic Index and % Genome Modified.

# Function that computes Spearman correlation on a resampled dataset
# (used as the `statistic` for boot() below).
spearman_boot_gi <- function(data, indices) {
  d <- data[indices, ]
  cor(d$Accuracy, d$gi, method = "spearman")
}

# Wrapper to run the bootstrap (BCa CI) for one technology/bin subset,
# correlating Accuracy against Genomic Index.
bootstrap_correlation_gi <- function(df, n_boot = 2000) {
  set.seed(42)  # for reproducibility

  boot_result <- boot(data = df, statistic = spearman_boot_gi, R = n_boot)

  ci <- boot.ci(boot_result, type = "bca")
  tibble(
    correlation = boot_result$t0,
    ci_lower = ci$bca[4],
    ci_upper = ci$bca[5]
  )
}

# Same pair as above, but correlating Accuracy against % Genome Modified
# ("val") instead of Genomic Index ("gi") - used in section 9b below.
spearman_boot_gm <- function(data, indices) {
  d <- data[indices, ]
  cor(d$Accuracy, d$val, method = "spearman")
}

bootstrap_correlation_gm <- function(df, n_boot = 2000) {
  set.seed(42)  # for reproducibility

  boot_result <- boot(data = df, statistic = spearman_boot_gm, R = n_boot)

  ci <- boot.ci(boot_result, type = "bca")
  tibble(
    correlation = boot_result$t0,
    ci_lower = ci$bca[4],
    ci_upper = ci$bca[5]
  )
}

#### 9a. Accuracy vs. Genomic Index ----
gi <- read.csv(paste0(getwd(), "/results/Genomic_Index/genomicIndexLMS.csv")) %>% 
  dplyr::select(c(Case, Bin, type, gi)) %>%
  dplyr::filter(type != "SNP")
acc <- read.csv(paste0(getwd(), "/results/Accuracy/accuracyLMS.csv")) %>% 
  dplyr::select(c(Cases, Accuracy, Bin_Size, Technology)) %>% 
  dplyr::rename(Case = Cases) %>%
  dplyr::rename(Bin = Bin_Size) %>%
  dplyr::rename(type = Technology)
acc[acc == "Sesame"] <- "SeSAMe"
data <- dplyr::inner_join(gi, acc, by = c("Case", "Bin", "type"))

correlations <- data %>%
  group_by(Bin, type) %>%
  group_modify(~ bootstrap_correlation_gi(.x)) %>%
  ungroup()
mat <- matrix(NA, nrow = 4, ncol = 3)
colnames(mat) <- unique(correlations$type)
rownames(mat) <- unique(correlations$Bin)
for (i in 1:nrow(correlations)) {
  row_index <- match(correlations$Bin[i], rownames(mat))
  col_index <- match(correlations$type[i], colnames(mat))
  mat[row_index, col_index] <- as.numeric(correlations$correlation[i])
}

colors <- rev(brewer.pal(n = 7, name = "RdBu"))
labels_matrix <- matrix(sprintf("%.3f", mat), 
                        nrow = nrow(mat), 
                        ncol = ncol(mat))
ph <- pheatmap(
  mat,
  color = colors,
  scale = "row", # Scales the values in each row/column/none to a z-score
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  display_numbers = labels_matrix,
  main = "Correlation Heatmap: Accuracy vs. Genomic Index",
  filename = paste0(getwd(), "/results/Correlations/Accuracy_vs_GenomicIndex_correlation.pdf"),
  width = 7,
  height = 5
)

# No worthy data from correlation for LMs.
#### 9b. Accuracy vs. Genome Modified ----
genome_modified <- read.csv(paste0(getwd(), "/results/Genome_modified/genome_modifiedLMS.csv")) %>% 
  dplyr::select(c(Case, Bin, Tech, val)) %>%
  dplyr::filter(Tech != "SNP") %>%
  dplyr::mutate(Bin = as.numeric(Bin))
acc <- read.csv(paste0(getwd(), "/results/Accuracy/accuracyLMS.csv")) %>% 
  dplyr::select(c(Cases, Accuracy, Bin_Size, Technology)) %>% 
  dplyr::rename(Case = Cases) %>%
  dplyr::rename(Bin = Bin_Size) %>%
  dplyr::rename(Tech = Technology)
data <- dplyr::inner_join(genome_modified, acc, by = c("Case", "Bin", "Tech"))

correlations <- data %>%
  group_by(Bin, Tech) %>%
  group_modify(~ bootstrap_correlation_gm(.x)) %>%
  ungroup()
mat <- matrix(NA, nrow = 4, ncol = 3)
colnames(mat) <- unique(correlations$Tech)
rownames(mat) <- unique(correlations$Bin)
for (i in 1:nrow(correlations)) {
  row_index <- match(correlations$Bin[i], rownames(mat))
  col_index <- match(correlations$Tech[i], colnames(mat))
  mat[row_index, col_index] <- as.numeric(correlations$correlation[i])
}

colors <- rev(brewer.pal(n = 7, name = "RdBu"))
labels_matrix <- matrix(sprintf("%.3f", mat), 
                        nrow = nrow(mat), 
                        ncol = ncol(mat))
ph <- pheatmap(
  mat,
  color = colors,
  scale = "row", # Scales the values in each row/column/none to a z-score
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  display_numbers = labels_matrix,
  main = "Correlation Heatmap: Accuracy vs. Genome Modified",
  filename = paste0(getwd(), "/results/Correlations/Accuracy_vs_GenomeModified_correlation.pdf"),
  width = 7,
  height = 5
)

#### 9c. SNP GI v. Methylation GI ----
gi <- read.csv(paste0(getwd(), "/results/Genomic_Index/genomicIndexLMS.csv")) %>% 
  dplyr::select(c(Case, Bin, type, gi)) %>%
  dplyr::filter(type != "SNP")
giSNP <- read.csv(paste0(getwd(), "/results/Genomic_Index/genomicIndexLMS.csv")) %>% 
  dplyr::select(c(Case, Bin, type, gi)) %>%
  dplyr::filter(type == "SNP") %>%
  dplyr::rename(giSNP = gi)
gi_combined <- inner_join(gi, giSNP, by = c("Case", "Bin")) %>% 
  dplyr::select(-c(type.y)) %>%
  dplyr::rename(Technology = type.x) %>%
  dplyr::rename(Genomic_Index = gi) %>%
  dplyr::rename(Genomic_Index_SNP = giSNP)
crr <- gi_combined %>% group_by(Technology) %>% 
  summarize(corr_coeff_pearson = cor(Genomic_Index, Genomic_Index_SNP, use = "complete.obs", method = "pearson"),
            corr_coeff_spearman = cor(Genomic_Index, Genomic_Index_SNP, use = "complete.obs", method = "spearman"),
            corr_coeff_kendall = cor(Genomic_Index, Genomic_Index_SNP, use = "complete.obs", method = "kendall"))
write.csv(crr, "results/Correlations/SNP_Methylation_GI_correlation.csv")

#### 9d. SNP GM v. Methylation GM ----
GenomeModifiedLMS <- read.csv(paste0(getwd(), "/results/Genome_modified/genome_modifiedLMS.csv")) %>%
  dplyr::select(-c("X")) %>%
  dplyr::filter(Tech != "SNP")

GenomeModifiedLMS_SNP <- read.csv(paste0(getwd(), "/results/Genome_modified/genome_modifiedLMS.csv")) %>% 
  dplyr::select(-c(X)) %>%
  dplyr::filter(Tech == "SNP") %>%
  dplyr::rename(valSNP = val) 

GenomeModified <- inner_join(GenomeModifiedLMS, GenomeModifiedLMS_SNP, by = c("Case")) %>%
  dplyr::select(-c(Tech.y, Bin.y)) %>%
  dplyr::rename(Technology = Tech.x, Bin = Bin.x) %>%
  group_by(Technology, Bin) %>%
  summarize(corr_coeff_pearson = cor(val, valSNP, use = "complete.obs", method = "pearson"),
            distrib_methyl = unname(shapiro.test(val)$statistic["W"]),
            distrib_SNP = unname(shapiro.test(valSNP)$statistic["W"]),
            p_val_SNP = unname(shapiro.test(valSNP)$p.value),
            p_val_methyl = unname(shapiro.test(val)$p.value))

write.csv(GenomeModified, paste0(getwd(), "/results/Correlations/SNP_Methylation_GenomeChanged_correlation.csv"))

#### 9e. Genomic Index correlation matrix (all callers incl. SNP) ----
# Pairwise Pearson r between each caller's 1 Mb Genomic Index and the
# SNP-array Genomic Index (ground truth), so it's visible at a glance
# which callers track SNP closely and which don't. LM is skipped here
# (see 9a's note) since most cases are missing a Genomic Index value
# for at least one caller, leaving too few complete pairs to correlate.
gi_wide <- read.csv(paste0(getwd(), "/results/Genomic_Index/genomicIndexLMS.csv")) %>%
  dplyr::select(Case, Bin, type, gi) %>%
  dplyr::filter(Bin == 1e+06) %>%
  tidyr::pivot_wider(id_cols = Case, names_from = type, values_from = gi)

gi_cor <- cor(dplyr::select(gi_wide, -Case), method = "pearson", use = "pairwise.complete.obs")

write.csv(gi_cor, paste0(getwd(), "/results/Correlations/GenomicIndex_correlation_matrix.csv"))

pdf(paste0(getwd(), "/results/Correlations/GenomicIndex_correlation_matrix.pdf"), width = 7, height = 7)
corrplot(
  gi_cor,
  method = "color",
  type = "upper",
  order = "original",
  addCoef.col = "black",
  number.digits = 3,
  tl.col = "black",
  tl.srt = 45,
  col = rev(brewer.pal(n = 7, name = "RdBu")),
  col.lim = c(-1, 1),
  title = "Genomic Index Correlation Matrix (1 Mb, Pearson r)",
  mar = c(0, 0, 2, 0)
)
dev.off()
