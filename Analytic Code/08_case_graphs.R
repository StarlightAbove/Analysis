### ============================================================
### 8. CASE GRAPHS
### ============================================================
# Per-case whole-genome CNV plots, saved as PDFs.
# Writes: <getwd()>/results/case_graphs/{LMS,LM,Normals}/...

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and results/), since all data paths
# are built from getwd().
source(paste0(getwd(), "/Analytic Code/00_setup.R"))
source(paste0(getwd(), "/Analytic Code/01_functions.R"))
source(paste0(getwd(), "/Analytic Code/02_load_data.R"))

# In the original single-file script, section 6 (Gene Concordance) reassigned
# `bins` to character values, and this section inherited that. Reproduced here
# so this script behaves exactly as before when run on its own (it affects the
# row order of the results and makes Bin_Size a discrete axis on the plots).
bins <- c("50000", "10000", "1e+05", "1e+06")

# Per-case, per-bin-size whole-genome CNV plots saved to disk as PDFs.
# All three cohorts get the plain plot_cnv_segments() track — LMS cases
# previously also got a gene-correlation heatmap side panel and gene-name
# notations on the log2 plot via geneAnno(), but that's been dropped.

#### 8a. LMS ----
for(c in LMS_cases){
  tag <- LMS_tags[[as.character(c)]]
  for(b in bins){
      my_plot <- plot_cnv_segments(df = rbind(labLMSProc(c, "MethylMaster", b),
                                              labLMSProc(c, "Sesame", b),
                                              labLMSProc(c, "Conumee", b)),
                                   title = paste0("CNV Segments — ", tag, " (LMS, ", b, " bp bins)"))
      dir.create(paste0(getwd(), "/results/case_graphs/LMS/",tag), showWarnings = FALSE, recursive = TRUE)
      ggsave(filename = paste0(getwd(), "/results/case_graphs/LMS/",tag,"/", tag,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)

  }
}

#### 8b. LM ----
for(c in LM_cases){
  tag <- LM_tags[[as.character(c)]]
  for(b in bins){
    my_plot <- plot_cnv_segments(df = rbind(LMStt(STT = c, tech = "MethylMaster", bin = b),
                                            LMStt(STT = c, tech = "Sesame", bin = b),
                                            LMStt(STT = c, tech = "Conumee", bin = b)),
                                 title = paste0("CNV Segments — ", tag, " (LM, ", b, " bp bins)"))
    dir.create(paste0(getwd(), "/results/case_graphs/LM/",tag), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = paste0(getwd(), "/results/case_graphs/LM/",tag,"/", tag,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)
  }
}

#### 8c. Normals ----
for(c in Normals){
  tag <- Normal_tags[[c]]
  for(b in bins){
    my_plot <- plot_cnv_segments(df = rbind(labNmrlProc(Sentrix = c, Technology = "MethylMaster", binSize = b),
                                            labNmrlProc(Sentrix = c, Technology = "Sesame", binSize = b),
                                            labNmrlProc(Sentrix = c, Technology = "Conumee", binSize = b)),
                                 title = paste0("CNV Segments — ", tag, " (Normal, ", b, " bp bins)"))
    dir.create(paste0(getwd(), "/results/case_graphs/Normals/",tag), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = paste0(getwd(), "/results/case_graphs/Normals/",tag,"/", tag,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)
  }
}
