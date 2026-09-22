### ============================================================
### 0. PACKAGE INSTALLATION
### ============================================================
# Run once on a new machine. Not sourced by any other script.

devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")

install.packages("BiocManager")

BiocManager::install(c(
  "GenomeInfoDb",
  "AnnotationDbi",
  "AnnotationHub",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "rentrez",
  "GenomicRanges",
  "GenomicFeatures",
  "pheatmap",
  "sesame"
))

install.packages(c(
  "tidyverse",
  "dplyr",
  "ggplot2",
  "patchwork",
  "cowplot",
  "gridGraphics",
  "gridExtra",
  "ggplotify",
  "readxl",
  "RColorBrewer",
  "reshape2",
  "purrr",
  "devtools",
  "corrplot"
))
