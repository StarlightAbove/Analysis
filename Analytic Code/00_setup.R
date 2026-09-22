### ============================================================
### 0. LIBRARY LOADING
### ============================================================
# Loads all packages. Sourced at the top of every analysis script.

# Bioconductor packages first...
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

# ...then tidyverse/dplyr last, so they win any function-name masking conflicts.
library(tidyverse)
library(dplyr)
library(ggplot2)
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
