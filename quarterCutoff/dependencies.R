devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
# Install BiocManager first to handle Bioconductor packages
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "GenomeInfoDb",
  "AnnotationDbi",
  "AnnotationHub",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "rentrez",
  "GenomicRanges",
  "GenomicFeatures",
  # "pheatmap",
  "sesame"
))

# Install CRAN packages
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
library(ggrepel)

# Then tidyverse and dplyr last so they win all masking conflicts
library(tidyverse)
library(dplyr)
library(ggplot2)
# library(ggrepel)
library(patchwork)
library(cowplot)
library(gridGraphics)
library(gridExtra)
library(ggplotify)
library(readxl)
library(RColorBrewer)
library(reshape2)
library(purrr)
library(corrplot)
library(boot)

