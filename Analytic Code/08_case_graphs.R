### ============================================================
### 8. CASE GRAPHS
### ============================================================
# Per-case whole-genome CNV plots, saved as PDFs.
# Writes: ~/Work/Analysis/quarterCutoff/case_graphs/{LMS,LM,Normals}/...

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

# Per-case, per-bin-size whole-genome CNV plots saved to disk as PDFs
# (LMS cases get the gene-annotated + heatmap version via geneAnno();
# LM cases and Normals get the plain plot_cnv_segments() track).

#### 8a. LMS ----
Genes <- c("MYC", "MYOCD", "CCNE1", "CDKN2A", "PTEN", "RB1", "TP53") 
for(c in LMS_cases){
  for(b in bins){
      my_plot <- geneAnno(Gene = Genes, 
                          db = rbind(labLMSProc(c, "MethylMaster", b), 
                                     labLMSProc(c, "Sesame", b), 
                                     labLMSProc(c, "Conumee", b)))
      ggsave(filename = paste0("~/Work/Analysis/quarterCutoff/case_graphs/LMS/",c,"/", c,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)
      
  }
}

#### 8b. LM ----
for(c in LM_cases){
  for(b in bins){
    my_plot <- plot_cnv_segments(df = rbind(LMStt(STT = c, tech = "MethylMaster", bin = b), 
                                            LMStt(STT = c, tech = "Sesame", bin = b), 
                                            LMStt(STT = c, tech = "Conumee", bin = b)))
    dir.create(paste0("~/Work/Analysis/quarterCutoff/case_graphs/LM/",c), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = paste0("~/Work/Analysis/quarterCutoff/case_graphs/LM/",c,"/", c,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)
  }
}

#### 8c. Normals ----
normals <- read.csv(
  "~/Work/Analysis/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")$Basename
for(c in normals){
  for(b in bins){
    my_plot <- plot_cnv_segments(df = rbind(labNmrlProc(Sentrix = c, Technology = "MethylMaster", binSize = b), 
                                            labNmrlProc(Sentrix = c, Technology = "Sesame", binSize = b),
                                            labNmrlProc(Sentrix = c, Technology = "Conumee", binSize = b)))
    dir.create(paste0("~/Work/Analysis/quarterCutoff/case_graphs/Normals/",c), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = paste0("~/Work/Analysis/quarterCutoff/case_graphs/Normals/",c,"/", c,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)
  }
}
