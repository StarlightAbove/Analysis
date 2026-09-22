### ============================================================
### 7. ACCURACY
### ============================================================
# Accuracy of each caller vs. SNP-array ground truth, plus CNV-only accuracy.
# Writes: quarterCutoff/Accuracy/accuracyLMS.csv, accuracyLM.csv, CNVaccuracyLMS.csv, CNVaccuracyLM.csv

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and quarterCutoff/), since all data paths
# are built from getwd().
source("00_setup.R")
source("01_functions.R")
source("02_load_data.R")

# In the original single-file script, section 6 (Gene Concordance) reassigned
# `bins` to character values, and this section inherited that. Reproduced here
# so this script behaves exactly as before when run on its own (it affects the
# row order of the results and makes Bin_Size a discrete axis on the plots).
bins <- c("50000", "10000", "1e+05", "1e+06")

# Per-chromosome true/false positive/negative accuracy of each caller
# against the SNP-array ground truth (see fpCheck()), plus an
# alternative "CNV-only" accuracy metric that ignores true negatives.

#### 7a. LMS ----
tech <- c("MethylMaster", "Sesame", "Conumee")
df <- expand.grid(
  Cases = LMS_cases,
  Bin_Size  = bins,
  Technology = tech,
  stringsAsFactors = FALSE
) %>% dplyr::mutate(Accuracy = 0)

for(c in LMS_cases){
  for(b in bins){
    for(t in tech){
      cs <- fpCheck(labLMSProc(c, t, b))
      cs <- cs[[2]]
      print(cs)
      idx <- which(df$Cases == c & 
                   df$Bin_Size  == b  & 
                   df$Technology == t)
      df$Accuracy[idx] <- cs
    }
  }
}

ggplot() +
  # Main lines for the three tools
  geom_line(
    data = df,
    aes(x = Bin_Size, y = Accuracy, colour = Technology, group = Technology),
    linewidth = 0.8
  ) +
  geom_point(
    data = df,
    aes(x = Bin_Size, y = Accuracy, colour = Technology),
    size = 2
  ) +
  # Facet by Case
  facet_wrap(~ Cases, ncol = 4) +
  # Scales and labels
  scale_colour_brewer(palette = "Set1", name = "Tool") +
  labs(
    title = "Accuracy across genome & cases",
    x     = "Bin size",
    y     = "Accuracy"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )
write.csv(df, "quarterCutoff/Accuracy/accuracyLMS.csv")

#### 7b. LM ----
tech <- c("MethylMaster", "Sesame", "Conumee")
df <- expand.grid(
  Cases = LM_cases,
  Bin_Size  = bins,
  Technology = tech,
  stringsAsFactors = FALSE
) %>% dplyr::mutate(Accuracy = 0)

for(c in LM_cases){
  for(b in bins){
    for(t in tech){
      cs <- fpCheck(LMStt(STT = c, bin = b, tech = t))
      cs <- cs[[2]]
      print(cs)
      idx <- which(df$Cases == c & 
                     df$Bin_Size  == b  & 
                     df$Technology == t)
      df$Accuracy[idx] <- cs
    }
  }
}

ggplot() +
  # Main lines for the three tools
  geom_line(
    data = df,
    aes(x = Bin_Size, y = Accuracy, colour = Technology, group = Technology),
    linewidth = 0.8
  ) +
  geom_point(
    data = df,
    aes(x = Bin_Size, y = Accuracy, colour = Technology),
    size = 2
  ) +
  # Facet by Case
  facet_wrap(~ Cases, ncol = 4) +
  # Scales and labels
  scale_colour_brewer(palette = "Set1", name = "Tool") +
  labs(
    title = "Accuracy across genome & cases",
    x     = "Bin size",
    y     = "Accuracy"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )
write.csv(df, "~/Work/Analysis/quarterCutoff/Accuracy/accuracyLM.csv")

#### 7c. Alternative LM accuracy calculation (CNV-only accuracy) ----
df <- expand.grid(
  Cases = LM_cases,
  Bin_Size  = bins,
  Technology = tech,
  stringsAsFactors = FALSE
) %>% dplyr::mutate(CNVAccuracy = 0)

for(c in LM_cases){
  for(b in bins){
    for(t in tech){
      cs <- fpCheck(LMStt(STT = c, bin = b, tech = t))
      acc <- sum(cs[[1]]$CNV_Only_Accuracy)/22
      idx <- which(df$Cases == c & 
                     df$Bin_Size  == b  & 
                     df$Technology == t)
      df$CNVAccuracy[idx] <- acc
    }
  }
}
write.csv(df, "quarterCutoff/Accuracy/CNVaccuracyLM.csv")

ggplot() +
  # Main lines for the three tools
  geom_line(
    data = df,
    aes(x = factor(Bin_Size), y = CNVAccuracy, colour = Technology, group = Technology),
    linewidth = 0.8
  ) +
  geom_point(
    data = df,
    aes(x = factor(Bin_Size), y = CNVAccuracy, colour = Technology),
    size = 2
  ) +
  # Facet by Case
  facet_wrap(~ Cases, ncol = 4) +
  # Scales and labels
  scale_colour_brewer(palette = "Set1", name = "Tool") +
  labs(
    title = "Accuracy across genome & cases",
    x     = "Bin size",
    y     = "Accuracy"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

#### 7d. Alternative LMS accuracy calculation (CNV-only accuracy) ----
df <- expand.grid(
  Cases = LMS_cases,
  Bin_Size  = bins,
  Technology = tech,
  stringsAsFactors = FALSE
) %>% dplyr::mutate(Accuracy = 0)

for(c in LMS_cases){
  for(b in bins){
    for(t in tech){
      cs <- fpCheck(labLMSProc(c, t, b))
      acc <- sum(cs[[1]]$CNV_Only_Accuracy)/22
      idx <- which(df$Cases == c & 
                     df$Bin_Size  == b  & 
                     df$Technology == t)
      df$CNVAccuracy[idx] <- acc
    }
  }
}

write.csv(df, "quarterCutoff/Accuracy/CNVaccuracyLMS.csv")

ggplot() +
  # Main lines for the three tools
  geom_line(
    data = df,
    aes(x = factor(Bin_Size), y = CNVAccuracy, colour = Technology, group = Technology),
    linewidth = 0.8
  ) +
  geom_point(
    data = df,
    aes(x = factor(Bin_Size), y = CNVAccuracy, colour = Technology),
    size = 2
  ) +
  # Facet by Case
  facet_wrap(~ Cases, ncol = 4) +
  # Scales and labels
  scale_colour_brewer(palette = "Set1", name = "Tool") +
  labs(
    title = "Accuracy across genome & cases",
    x     = "Bin size",
    y     = "Accuracy"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )
