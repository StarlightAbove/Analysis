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
  "pheatmap",
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
  "devtools"
))
install.packages("corrplot")


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

## Functions ----
labLMSProc <- function(STTq, Technology, binSize){
  
  
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/correlative.csv")) %>% filter(STT == STTq)
  
  cnvMatch <- read_delim(paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_LMS_Probe_and_segment_level_data_01Feb2026/STT", 
                                STTq, "_Recentered_Segment_level_data_01Feb2026.segment.txt"))
  
  # methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
  #paste0(correlationSheet$Sentrix_ID[1],"_", 
  #     correlationSheet$Sentrix_Position[1], "/"), 
  # "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0(getwd(), "/Outputs/Conumee/LabLMS/", binSize,"/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                                 correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0(getwd(), "/Outputs/SeSAMe/LabLMS/", binSize,"/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
                                                                                               correlationSheet$Sentrix_Position[1],".csv")))
  
  methylMatch <- read.csv(paste0(getwd(), "/Outputs/MethylMaster/LabLMS/", binSize,"/", 
                                 correlationSheet$Sentrix_ID[1],"_", correlationSheet$Sentrix_Position[1],"/autocorrected_regions.csv"))
  
  cnvMatch <- cnvMatch %>% dplyr::filter(!(Type == "LOH")) %>% dplyr::filter(Chromosome != 24 & Chromosome != 25) %>% dplyr::select("Chromosome", "StartPosition", "StopPosition", "Median Log2 Ratio") %>% 
    dplyr::rename(chrom = "Chromosome", loc.start = "StartPosition", loc.end = "StopPosition", seg.mean = "Median Log2 Ratio") %>%
    dplyr::mutate(seg.mean = as.numeric(seg.mean)) %>%
    dplyr::mutate(CNVStatus = case_when(
      seg.mean <= -0.25 ~ "Deletion", 
      seg.mean >= 0.25 ~ "Amplification", 
      TRUE ~ "Normal"
    ), type = "SNP")
  
  
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.25 ~ "Deletion", 
        seg.mean >= 0.25 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- rbind(methylMatch, cnvMatch) %>% mutate(Gene = as.character(0)) %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.25 ~ "Deletion", 
      seg.mean >= 0.25 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- rbind(conumeeMatch, cnvMatch) %>% mutate(Gene = as.character(0)) %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.25 ~ "Amplification",
                                                                                 seg.mean < -0.25 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- rbind(sesameMatch, cnvMatch) %>% mutate(Gene = as.character(0), seg.mean = as.numeric(seg.mean)) %>% arrange(chrom)
  }
  combinedSet
}
LMStt <- function(STT, bin, tech){
  
  corrSheet <- read.csv(paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")) %>%
    dplyr::filter(STT == !!STT)
  sentrixID <- corrSheet$Sentrix_ID[1]
  sentrixPos <- corrSheet$Sentrix_Position[1]
  
  SNP <- read_delim(paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_LM_Probe_and_segment_level_data_01Feb2026/STT",
                           STT,"_Segment_level_data_01Feb2026.segment.txt"))
  SNP <- SNP %>% dplyr::filter(!(Type == "LOH")) %>% dplyr::filter(Chromosome != 24 & Chromosome != 25) %>% dplyr::select("Chromosome", "StartPosition", "StopPosition", "Median Log2 Ratio") %>% 
    dplyr::rename(chrom = "Chromosome", loc.start = "StartPosition", loc.end = "StopPosition", seg.mean = "Median Log2 Ratio") %>%
    dplyr::mutate(seg.mean = as.numeric(seg.mean)) %>%
    dplyr::mutate(CNVStatus = case_when(
      seg.mean <= -0.25 ~ "Deletion", 
      seg.mean >= 0.25 ~ "Amplification", 
      TRUE ~ "Normal"
    ), type = "SNP")
  
  
  # MethylMasteR
  cnvMethyl <- NULL;
  if(tech == "MethylMaster"){
    outputDirMethyl <- paste0(getwd(), "/Outputs/MethylMaster/LM/",bin,"/",sentrixID,"_",sentrixPos,"/autocorrected_regions.csv")
    cnvMethyl <- read.csv(outputDirMethyl) %>%
      dplyr::select(c("Chromosome", "bp.Start", "bp.End", "Mean")) %>% dplyr::rename(
        chrom = "Chromosome",
        loc.start = "bp.Start",
        loc.end = "bp.End",
        seg.mean = "Mean"
      ) %>% dplyr::mutate(
        CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                              seg.mean < -0.3 ~ "Deletion",
                              TRUE ~ "Normal"), type = "MethylMaster"
      )
    cnvMethyl$chrom <- as.numeric(str_remove_all(cnvMethyl$chrom, pattern = "chr"))
    cnvMethyl <- cnvMethyl %>% filter(!(is.na(chrom))) %>% arrange(chrom)
  }
  
  
  # SeSAMe
  sesameOutput <- NULL
  if(tech == "Sesame"){
    outputDirSesame <- paste0(getwd(),"/Outputs/SeSAMe/LM/bins/",bin,"/","segments_",sentrixID, "_", sentrixPos, ".csv")
    sesameOutput <- read.csv(outputDirSesame) %>% dplyr::select(c("chrom", "loc.start", 
                                                                  "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                                    chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                                      !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                                        chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                                          CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                                seg.mean < -0.3 ~ "Deletion",
                                                                                                TRUE ~ "Normal"), type = "SeSAMe") 
  }
  
  
  # Conumee
  case <- NULL
  if(tech == "Conumee"){
    outputDirConumee <- paste0(getwd(), "/Outputs/Conumee/LMData/bins/",bin,"/",sentrixID,"_",sentrixPos,".csv")
    case <- read.csv(outputDirConumee)
    case <- case %>% dplyr::select(-c("ID", "bstat")) %>% 
      mutate(CNVStatus = case_when(case$seg.mean > 0.25 ~ "Amplification",
                                   case$seg.mean < -0.25 ~ "Deletion",
                                   TRUE ~ "Normal"), 
             chrom = as.numeric(str_remove_all(chrom, "chr"))) %>% 
      dplyr::arrange(chrom) %>% mutate(type = "Conumee") %>%
      dplyr::select(-c("num.mark", "pval", "seg.median", "X"))
  }
  
  output <- rbind(SNP,case,sesameOutput,cnvMethyl) %>% dplyr::mutate(Gene = as.character(0))
  output$chrom <- as.character(output$chrom)
  output
}
labNmrlProc <- function(Sentrix, Technology, binSize){
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(
    "/LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv") %>% filter(Basename == Sentrix)
  
  methylMatch <- read.csv(paste0(getwd(), "/Outputs/MethylMaster/Normals/", binSize, "/", Sentrix, "/autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0(getwd(),"/Outputs/Conumee/LabNormals/", binSize,"/", paste0(Sentrix,".csv")))
  
  sesameMatch <- read.csv(paste0(getwd(),"/Outputs/SeSAMe/Normals/", binSize,"/", paste0("segments_",Sentrix,".csv")))
  
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.25 ~ "Deletion", 
        seg.mean >= 0.25 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- methylMatch %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.25 ~ "Deletion", 
      seg.mean >= 0.25 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- conumeeMatch %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                 seg.mean < -0.3 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- sesameMatch %>% arrange(chrom)
  }
  
  
  
  combinedSet <- combinedSet %>% dplyr::mutate(Gene = as.character(0))
  combinedSet
}

caseCorr <- function(IDs, bin){
  
  Methyl <- c()
  Ses <- c()
  Con <- c()
  
  for(ID in IDs){
    MethylMaster <- read.csv(paste0(getwd(),"/Outputs/MethylMaster/Normals/",
                                    bin,"/",ID,"/autocorrected_regions.csv"))
    MethylMaster <- median(MethylMaster$Mean)
    Sesame <- read.csv(paste0(getwd(),"/Outputs/SeSAMe/Normals/",
                              bin,"/segments_",ID,".csv"))
    Sesame <- median(Sesame$seg.mean)
    Conumee <- read.csv(paste0(getwd(),"/Outputs/Conumee/LabNormals/",
                               bin,"/",ID,".csv"))
    Conumee <- median(Conumee$seg.mean)
    
    Methyl <- c(Methyl, MethylMaster)
    Ses <- c(Ses, Sesame)
    Con <- c(Con,Conumee)
  }
  
  print(Methyl)
  print(Ses)
  print(Con)
  
  acc <- data.frame(Conumee_Median = median(Con), 
                    Conumee_SD = sd(Con), 
                    Conumee_min = min(Con),
                    Conumee_max = max(Con),
                    Sesame_Median = median(Ses), 
                    Sesame_SD = sd(Ses), 
                    Sesame_min = min(Ses),
                    Sesame_max = max(Ses),
                    Methyl_Median = median(Methyl), 
                    Methyl_SD = sd(Methyl),
                    Methyl_min = min(Methyl),
                    Methyl_max = max(Methyl))
  acc
  
}
GenomicIndex <- function(dfCH3, stt = 0, bin, sw, LMSorLM = "LMS"){
  recentered <- NULL
  rc <- NULL
  if(LMSorLM == "LM"){
    recentered <- ""
    rc <- ""
  } else { 
    recentered <- "Recentered_"
    rc <- ".RC"
  }
  if(sw == F && stt != 0){
    dfSNP <- read_delim(paste0(getwd(), "/LabData/",LMSorLM,"_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_",LMSorLM,"_CNV_calls_01Feb2026/",recentered,"STT", stt, rc,".OSCHP.segments.txt")) %>%
      dplyr::filter(Chromosome != "X", `Size (kbp)` >= bin/1000) 
    print(dfSNP)
    modChrom <- length(unique(dfSNP$Chromosome))
    CNVCount <- nrow(dfSNP)
    ret <- (CNVCount^2)/modChrom
    return(ret)
  } else {
    # Filter out for clinically significant segments %>% filter(CNVStatus != Normal)
    dfCH3f <- dfCH3 %>% 
      dplyr::filter(CNVStatus != "Normal", chrom != "X", type != "SNP")
    # Count
    CNVCount <- nrow(dfCH3f)
    modChrom <- length(unique(dfCH3f$chrom))
    
    
    # Return value
    ret <- (CNVCount^2)/modChrom
    return(ret) 
  }
}

GenomicIndexIntermediateMatrix <- function(dfCH3, stt = 0, bin, LMSorLM = "LMS", tech){
  print(stt)
  recentered <- NULL
  rc <- NULL
  if(LMSorLM == "LM"){
    recentered <- ""
    rc <- ""
  } else { 
    recentered <- "Recentered_"
    rc <- ".RC"
  }
  dfSNP <- read_delim(paste0(getwd(), "/LabData/",LMSorLM,"_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_",LMSorLM,"_CNV_calls_01Feb2026/",recentered,"STT", stt, rc,".OSCHP.segments.txt")) %>%
    dplyr::filter(Chromosome != "X", `Size (kbp)` >= 1e+04) %>% # 10Mb cutoff seems to generate reasonable values
    dplyr::select(c("Median Log2Ratio", "Cytoband Start", "Cytoband End", 
                    "Chromosome", "Type")) %>%
    dplyr::rename(seg.mean = "Median Log2Ratio", loc.start = "Cytoband Start",
                  loc.end = "Cytoband End", chrom = Chromosome, CNVStatus = Type) %>%
    dplyr::mutate(type = "SNP", Gene = 0) %>%
    dplyr::filter(seg.mean > 0.25 | seg.mean < -0.25)

  
  dfCH3f <- dfCH3 %>% 
    dplyr::filter(CNVStatus != "Normal", chrom != "X", type != "SNP")

  output <- rbind(dfCH3f, dfSNP) %>% dplyr::mutate(Case = stt, Bin = bin)
  
  
  output
}

GenomeModified <- function(dfCH3){
  df <- dfCH3 %>% dplyr::filter(type != "SNP", CNVStatus != "Normal") %>%
    dplyr::mutate(width = loc.end - loc.start) 
  modifiedWidth <- sum(df$width)
  hg19_info <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(assembled == "TRUE")
  hg19_info <- hg19_info[1:22,]
  hg19_total <- sum(hg19_info$size)
  
  return((modifiedWidth/hg19_total) * 100)
  
}

fpCheck <- function(df) {
  library(GenomicRanges)
  library(dplyr)
  # Removes genes, it's unnecessary for this purpose
  
  pred_df <- df %>% filter(!(type == "SNP")) # The methylation dataset
  pred_df$chrom <- as.character(pred_df$chrom)
  truth_df <- df %>% filter(type == "SNP") # The SNP dataset
  truth_df$chrom <- as.character(truth_df$chrom)
  
  # Both chromosomes are made into characters because we want to use them as factors
  # ---- Step 1: hg19 chromosome sizes (numeric chromosomes 1 to 22) ----
  hg19_info <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(assembled == "TRUE") # Just gets the main dataset from UCSC.
  hg19_info <- hg19_info[1:22,] # Only uses the first 22 chromosomes.
  hg19_chr_sizes <- data.frame(
    Chromosome = as.character(1:22),
    Genome_bp = hg19_info$size
  ) # Makes them into a referenceable dataframe.
  
  # ---- Step 2: Convert to GRanges ---- Does what it says on the tin
  truth_gr <- GRanges(seqnames = truth_df$chrom,
                      ranges = IRanges(start = truth_df$loc.start, end = truth_df$loc.end),
                      CNV = truth_df$CNVStatus)
  
  pred_gr <- GRanges(seqnames = pred_df$chrom,
                     ranges = IRanges(start = pred_df$loc.start, end = pred_df$loc.end),
                     CNV = pred_df$CNVStatus)
  
  # ---- Step 3: Overlaps and TP ----
  hits <- findOverlaps(pred_gr, truth_gr) # Finds the overlaps between the two GRanges
  
  overlap_ranges <- pintersect(pred_gr[queryHits(hits)], truth_gr[subjectHits(hits)])
  pred_cnv <- mcols(pred_gr)$CNV[queryHits(hits)] # Returns whether the intersections are deletions, amplifications or normals
  truth_cnv <- mcols(truth_gr)$CNV[subjectHits(hits)] # Returns whether the intersections are deletions, amplifications or normals
  
  # The returns are for the corresponding datasets.
  
  tp_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = pred_cnv,
    truth_cnv = truth_cnv
  ) # Creates a dataframe from the given data
  tp_df <- tp_df %>%
    filter(pred_cnv == truth_cnv) %>%
    group_by(Chromosome) %>%
    summarise(TP_bp = sum(width), .groups = "drop")
  
  # ---- Step 4: False Positives ----
  fp_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = pred_cnv,
    truth_cnv = truth_cnv
  ) %>% filter(pred_cnv != truth_cnv, truth_cnv == "Normal") %>%
    group_by(Chromosome) %>%
    summarise(FP_bp = sum(width), .groups = "drop") # Calculates the false positives
  
  # ---- Step 5: False Negatives ----
  fn_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = pred_cnv,
    truth_cnv = truth_cnv
  ) %>% filter(pred_cnv != truth_cnv, pred_cnv == "Normal") %>%
    group_by(Chromosome) %>%
    summarise(FN_bp = sum(width), .groups = "drop") # Calculates the false negatives
  
  # ---- Step 6: Merge and compute accuracy ----
  final_eval <- hg19_chr_sizes %>%
    left_join(tp_df, by = "Chromosome") %>%
    left_join(fp_df, by = "Chromosome") %>%
    left_join(fn_df, by = "Chromosome") %>%
    mutate(across(c(TP_bp, FP_bp, FN_bp), ~replace_na(., 0))) %>%
    mutate(
      CNV_bp = TP_bp + FP_bp + FN_bp,
      TN_bp = Genome_bp - CNV_bp,
      Accuracy = (TP_bp + TN_bp) / Genome_bp,
      CNV_Only_Accuracy = ifelse((TP_bp + FP_bp + FN_bp) == 0, NA,
                                 TP_bp / (TP_bp + FP_bp + FN_bp))
    ) %>%
    arrange(as.numeric(Chromosome)) # Joins all the analysis together, so it can be seen at once, and then calculates the accuracy across each case
  
  accuracy <- median(final_eval$Accuracy)
  
  return(list(final_eval, accuracy))
}
NoGraphGeneGen <- function(Gene, db = NULL, case){ # Updated function.
  
  
  # This section annotates the processed dataframe with genes by
  # Querying Entrez to get the genes we want to, get their locations
  # Overlap those locations onto the existing dataframe, and voila
  # We know where the genes are!
  entrez_idsDB <- mapIds(org.Hs.eg.db, 
                         keys=Gene, 
                         column="ENTREZID", 
                         keytype="SYMBOL",
                         multiVals="first")
  chr_locations <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = entrez_idsDB,
                                         columns = c("CHR", "CHRLOC", "CHRLOCEND"),
                                         keytype = "ENTREZID")
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_ranges <- genes(txdb)
  target_genes_ranges <- gene_ranges[gene_ranges$gene_id %in% entrez_idsDB]
  coords_df <- as.data.frame(target_genes_ranges)
  coords_df$entrez_id <- target_genes_ranges$gene_id
  coords_final <- coords_df[, c("entrez_id", "seqnames", "start", "end", 
                                "strand", "width")] %>% 
    dplyr::mutate(chrom=seqnames) %>% 
    dplyr::mutate(chrom = as.numeric(gsub("chr","",chrom)), 
                  Gene = Gene) %>% 
    dplyr::select(-c(seqnames))
  rownames(coords_final) <- Gene
  
  
  
  
  if(!is.null(db)){
    print("is here")
    dt <- db %>% dplyr::filter((Gene == "0")) %>% 
      dplyr::mutate(Gene = NA) # We make sure we don't hold onto anything that
    # is not associated with a gene
    
    gcoords <- coords_final %>% 
      dplyr::mutate(loc.start = as.numeric(start), loc.end = as.numeric(end)) %>% 
      dplyr::mutate(seg.mean = 0, CNVStatus = "Normal", type = "Gene") %>%
      dplyr::mutate(Gene = Gene) %>% 
      dplyr::select(-c("strand", "width", "entrez_id", "start", "end"))
    
    grC <- makeGRangesFromDataFrame(gcoords, 
                                    start.field = "loc.start",
                                    end.field = "loc.end",
                                    keep.extra.columns = TRUE)
    grDt <- makeGRangesFromDataFrame(dt,
                                     start.field = "loc.start",
                                     end.field = "loc.end",
                                     keep.extra.columns = TRUE)
    
    overlaps <- findOverlaps(grDt, grC)
    
    # Extract intersection and bind gene annotation immediately before any row changes
    intersection_df <- as.data.frame(pintersect(
      grDt[queryHits(overlaps)], 
      grC[subjectHits(overlaps)]
    )) %>%
      dplyr::select(-c("strand", "width", "hit", "Gene")) %>%
      dplyr::rename(loc.start = start, loc.end = end) %>%
      dplyr::mutate(chrom = as.numeric(seqnames)) %>%
      dplyr::select(-c("seqnames"))
    
    gene_annotations <- coords_final[subjectHits(overlaps), "Gene", drop = FALSE]
    overlapping_regions <- cbind(intersection_df, Gene = gene_annotations$Gene)
    
    # Find genes with no overlap for each type
    all_types <- unique(dt$type)
    missing_rows <- list()
    
    for(t in all_types){
      # Detects if the gene is split between two or more segments by taking the average instead
      dt_type <- dt %>% dplyr::filter(type == t)
      grDt_type <- makeGRangesFromDataFrame(dt_type,
                                            start.field = "loc.start",
                                            end.field = "loc.end",
                                            keep.extra.columns = TRUE)
      matched_genes <- unique(overlapping_regions$Gene[overlapping_regions$type == t])
      missing_genes <- setdiff(rownames(coords_final), matched_genes)
      
      for(g in missing_genes){
        gene_row <- coords_final[g, ]
        gene_mid <- (gene_row$start + gene_row$end) / 2
        
        # Find flanking segments on same chromosome
        dt_chr <- dt_type %>% dplyr::filter(chrom == gene_row$chrom)
        
        if(nrow(dt_chr) == 0) next
        
        # Left flank: segments ending before gene midpoint
        left <- dt_chr %>% dplyr::filter(loc.end < gene_mid) %>% 
          dplyr::slice_max(loc.end, n = 1)
        # Right flank: segments starting after gene midpoint  
        right <- dt_chr %>% dplyr::filter(loc.start > gene_mid) %>% 
          dplyr::slice_min(loc.start, n = 1)
        
        if(nrow(left) == 0 || nrow(right) == 0) next
        
        avg_seg <- mean(c(left$seg.mean, right$seg.mean))
        avg_cnv <- ifelse(left$CNVStatus == right$CNVStatus, 
                          left$CNVStatus, "Normal")
        
        missing_rows[[length(missing_rows) + 1]] <- data.frame(
          chrom = gene_row$chrom,
          loc.start = gene_row$start,
          loc.end = gene_row$end,
          seg.mean = avg_seg,
          type = t,
          CNVStatus = avg_cnv,
          Gene = g,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if(length(missing_rows) > 0){
      overlapping_regions <- rbind(overlapping_regions, 
                                   dplyr::bind_rows(missing_rows))
    } # Makes sure we have gotten everything
    
    
    overlapping_regions <- unique(overlapping_regions[c(6, 1, 2, 3, 4, 5, 7)])
    overlapping_regions$type <- ifelse(
      overlapping_regions$type == "Conumee", "Gene_Conumee",
      ifelse(overlapping_regions$type == "SeSAMe", "Gene_SeSAMe", 
             ifelse(overlapping_regions$type == "MethylMaster", "Gene_MMasteR",
                    ifelse(overlapping_regions$type == "SNP", "Gene_SNP", "Unknown")))
    )
    overlapping_regions <- overlapping_regions %>% dplyr::mutate(case = case)
    overlapping_regions
    
  } else {
    coords_final
  }
}

plot_cnv_segments <- function(df) {
  
  df <- df %>%
    mutate(
      type_group = ifelse(type == "SNP", "SNP", "non-SNP")
    ) %>%
    arrange(chrom, loc.start)
  
  df$seg.mean <- as.numeric(df$seg.mean)
  
  # Calculate chromosome cumulative positions
  chr_lengths <- df %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarize(chr_len = max(loc.end), .groups = "drop") %>%
    arrange(chrom) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
    mutate(chr_mid = chr_start + chr_len / 2)
  
  df <- df %>%
    left_join(chr_lengths, by = "chrom") %>%
    mutate(
      start_cum = loc.start + chr_start,
      end_cum = loc.end + chr_start
    ) %>%mutate(type = ifelse(type == "SNP", "SNP", type)) %>% 
    mutate(type = as.factor(type))
  
  # Vertical chromosome boundaries
  chr_boundaries <- chr_lengths %>%
    mutate(x = chr_start) %>%
    dplyr::select(chrom, x)
  
  x_breaks <- chr_lengths$chr_mid
  x_labels <- paste0("chr", chr_lengths$chrom)
  x_labels <- c(x_labels)
  print(x_labels)
  print(df)
  df$Gene[df$Gene == "0"] <- NA
  gene_labels <- subset(df, type == "Gene_SNP" & !is.na(Gene) & Gene != "0") %>%
    arrange(start_cum) %>%
    mutate(
      too_close = c(FALSE, diff(start_cum) < 5e7),
      y_offset  = ifelse(too_close, -Inf, -Inf)  # both bottom
    )
  # Plot
  p <- ggplot(df, aes(x = start_cum, xend = end_cum, y = seg.mean, yend = seg.mean, label = Gene)) +
    geom_segment(aes(color = type), size = 0.7, alpha = 0.8) +
    geom_point(data = subset(df, type == "Gene_Conumee"), aes(x = start_cum, y = seg.mean), color = "green", size = 3) + 
    geom_point(data = subset(df, type == "Gene_SNP"), aes(x = start_cum, y = seg.mean), color = "gray40", size = 3) + 
    geom_point(data = subset(df, type == "Gene_SeSAMe"), aes(x = start_cum, y = seg.mean), color = "red", size = 3) + 
    geom_point(data = subset(df, type == "Gene_MMasteR"), aes(x = start_cum, y = seg.mean), color = "blue", size = 3) + 
    geom_vline(data = gene_labels,
               aes(xintercept = start_cum),
               color = "orange", linetype = "solid", 
               linewidth = 0.4, alpha = 0.6) +
    geom_label(data = gene_labels,
               aes(x = start_cum, y = -Inf, label = Gene),
               vjust  = ifelse(gene_labels$too_close, 1, 0),
               hjust = 0,
               angle = 90,
               size = 3,
               fill = "white",
               color = "orange",
               label.size = 0.2,
               inherit.aes = FALSE) +
    geom_vline(data = chr_boundaries, aes(xintercept = x), color = "grey70", linetype = "dashed") +
    # scale_color_manual(values = c("Amplification" = "red", "Deletion" = "blue", "Normal" = "black")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("SeSAMe" = "red", "MethylMaster" = "blue", 
                                  "Conumee" = "green", "SNP" = "black", "Gene" = "orange")) +
    labs(
      x = "Genomic Position (across chromosomes)",
      y = "Segment Mean (log2 ratio)",
      title = "CNV Segments Across Genome"
    ) + geom_hline(yintercept = -0.25, linetype = "dotted", color = "black") + 
    geom_hline(yintercept = 0.25, linetype = "dotted", color = "black") + 
    scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0, 0.25,0.5,0.75)) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(5, 5, 40, 5, "pt")
    )
  
  return(p)
}
geneGen <- function(Gene, db = NULL){
  entrez_idsDB <- mapIds(org.Hs.eg.db, 
                         keys=Gene, 
                         column="ENTREZID", 
                         keytype="SYMBOL",
                         multiVals="first")
  chr_locations <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = entrez_idsDB,
                                         columns = c("CHR", "CHRLOC", "CHRLOCEND"),
                                         keytype = "ENTREZID")
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  gene_ranges <- genes(txdb)
  target_genes_ranges <- gene_ranges[gene_ranges$gene_id %in% entrez_idsDB]
  coords_df <- as.data.frame(target_genes_ranges)
  coords_df$entrez_id <- target_genes_ranges$gene_id
  coords_final <- coords_df[, c("entrez_id", "seqnames", "start", "end", 
                                "strand", "width")] %>% 
    dplyr::mutate(chrom=seqnames) %>% 
    dplyr::mutate(chrom = as.numeric(gsub("chr","",chrom)), 
                  Gene = Gene) %>% 
    dplyr::select(-c(seqnames))

  rownames(coords_final) <- Gene
  
  if(!is.null(db)){
    print("is here")
    dt <- db %>% dplyr::filter((Gene == "0")) %>% 
      dplyr::mutate(Gene = NA)
    gcoords <- coords_final %>% 
      dplyr::mutate(loc.start = as.numeric(start), loc.end = as.numeric(end)) %>% 
      dplyr::mutate(seg.mean = 0, CNVStatus = "Normal", type = "Gene") %>%
      dplyr::mutate(Gene = Gene) %>% 
      dplyr::select(-c("strand", "width", "entrez_id", "start", "end"))
    grC <- makeGRangesFromDataFrame(gcoords, 
                                    start.field = "loc.start",
                                    end.field = "loc.end",
                                    keep.extra.columns = TRUE)
    grDt <- makeGRangesFromDataFrame(dt,
                                     start.field = "loc.start",
                                     end.field = "loc.end",
                                     keep.extra.columns = TRUE)
    overlaps <- findOverlaps(grDt, grC)
    overlapping_regions <- pintersect(grDt[queryHits(overlaps)], grC[subjectHits(overlaps)])
    overlapping_regions <- as.data.frame(overlapping_regions) %>% 
      dplyr::select(-c("strand", "width", "hit", "Gene")) %>%
      dplyr::rename(loc.start = start, loc.end = end) %>%
      dplyr::mutate(chrom = as.numeric(seqnames)) %>%
      dplyr::select(-c("seqnames"))
    overlapping_regions <- dplyr::left_join(overlapping_regions, coords_final, 
                                            by=c("chrom"="chrom",
                                                 "loc.start" = "start", 
                                                 "loc.end" = "end")) %>%
      dplyr::select(-c("entrez_id", "strand", "width")) 
    overlapping_regions <- unique(overlapping_regions[c(6, 1, 2, 3, 4, 5, 7)])
    overlapping_regions$type <- ifelse(
      overlapping_regions$type == "Conumee", "Gene_Conumee",
      ifelse(overlapping_regions$type == "SeSAMe", "Gene_SeSAMe", 
             ifelse(overlapping_regions$type == "MethylMaster", "Gene_MMasteR",
                    ifelse(overlapping_regions$type == "SNP", "Gene_SNP", "Unknown")))
    )
    overlapping_regions <- overlapping_regions %>% drop_na()
    
    
    # Developing visualization
    mat <- matrix(NA, nrow = 4, ncol = length(unique(overlapping_regions$Gene)))
    overlapping_regions <- overlapping_regions %>% drop_na()
    colnames(mat) <- unique(overlapping_regions$Gene)
    rownames(mat) <- unique(overlapping_regions$type)
    for (i in 1:nrow(overlapping_regions)) {
      row_index <- match(overlapping_regions$type[i], rownames(mat))
      col_index <- match(overlapping_regions$Gene[i], colnames(mat))
      mat[row_index, col_index] <- as.numeric(overlapping_regions$seg.mean[i])
    }
    mat[is.na(mat)] <- 0
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
      show_rownames = TRUE, # Set to TRUE if you have few genes
      display_numbers = labels_matrix,
      main = "Gene Expression Heatmap (Log2 Ratio / Z-score)"
    )
    list(ph, overlapping_regions)
    
  } else {
    coords_final
  }
}
geneAnno <- function(Gene, db = NULL){
  reference <- geneGen(Gene = Gene, db = db)
  db <- rbind(db, reference[[2]])
  heatmap <- reference[[1]]
  pheatmap_ggplot <- as.ggplot(heatmap$gtable)
  plot <- plot_cnv_segments(db)
  return(plot + pheatmap_ggplot + plot_layout(widths = unit(c(20,8), c("null","null"))))
}


#### ------ READ ALL DATA ------ ####
# Getting the list of cases
LMS_cases <- read_csv(paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/correlative.csv"))
LMS_cases <- LMS_cases$STT

LM_cases <- read_csv(paste0(getwd(),"/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv"))
LM_cases <- LM_cases$STT

Normals <- read_csv("LabData/Normal_smooth_muscle_EPIC_data/idat_files/Sample_Sheet_Normal.csv")
Normals <- Normals$Basename

bins <- c(10000, 50000, 1e+05, 1e+06)

### Prelim Data ----
hg19_info <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(assembled == "TRUE")
hg19_info <- hg19_info[1:22,]
hg19_total <- sum(hg19_info$size)



### NORMALS ----

normalsFrame10kb <- caseCorr(IDs = Normals, bin = 10000) %>% dplyr::mutate(bin = 10000)
normalsFrame100kb <- caseCorr(IDs = Normals, bin = 1e+05) %>% dplyr::mutate(bin = 1e+05)
normalsFrameDef <- caseCorr(IDs = Normals, bin = 50000) %>% dplyr::mutate(bin = 50000)
normalsFrame1Mb <- caseCorr(IDs = Normals, bin = 1e+06) %>% dplyr::mutate(bin = 1e+06)

NormalsFrame <- rbind(normalsFrame10kb, normalsFrame100kb, normalsFrameDef, 
                      normalsFrame1Mb) %>% arrange(bin) 
write.csv(NormalsFrame, file = paste0(getwd(), "/quarterCutoff/normals.csv"))

# We can observe none of the bins have anything outside of the cutoff


### GENOMIC INDEX ----
## LMS ----
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

write.csv(summ, paste0(getwd(), "/quarterCutoff/Genomic_Index/genomicIndexLMS.csv"))

## LM ----
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

write.csv(LMGI, paste0(getwd(), "/quarterCutoff/Genomic_Index/genomicIndexLM.csv"))

## Plotting ----
giplot <- function(df){
  p <- ggplot(df, aes(x = Case, y = gi, fill = type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~ Bin, ncol = 2, labeller = labeller(Bin = c(
      `5e+04` = "Default",
      `1e+04` = "10 kb",
      `1e+05` = "100 kb",
      `1e+06` = "1 Mb"
    ))) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Genomic Index by Case and Tool",
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

df <- read.csv("~/Work/Analysis/quarterCutoff/Genomic_Index/genomicIndexLM.csv") %>%
  select(Case, Bin, type, gi) %>%
  mutate(Case = factor(Case)) %>%
  mutate(Bin = factor(Bin, levels = c(1e+04, 5e+04, 1e+05, 1e+06)))
LMPlot <- giplot(df)

df <- read.csv(paste0(getwd(), "/quarterCutoff/Genomic_Index/genomicIndexLMS.csv")) %>%
  select(Case, Bin, type, gi) %>%
  mutate(Case = factor(Case)) %>%
  mutate(Bin = factor(Bin, levels = c(1e+04, 5e+04, 1e+05, 1e+06)))
LMSPlot <- giplot(df)


### GENOME MODIFIED ----
## LMS ----
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
write.csv(geneMod, file = paste0(getwd(), "/quarterCutoff/Genome_modified/genome_modifiedLMS.csv"))
df <- read.csv("quarterCutoff/Genome_modified/genome_modifiedLMS.csv") %>% dplyr::select(-c("X"))
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

## LM ----
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
write.csv(LMdf, file = paste0(getwd(), "/quarterCutoff/Genome_modified/genome_modifiedLM.csv"))





### GENE CONCORDANCE ----
## LMS ----
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


# Saving raw data
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

# Looking for concordance
genes <- names(aps)
bins <- c("50000", "10000", "1e+05", "1e+06")
states <- c("Conumee–default bin size",
            "Conumee–10kb",
            "Conumee–100kb",
            "Conumee–1Mb",            
            "Sesame–default bin size",
            "Sesame–10kb",
            "Sesame–100kb",
            "Sesame–1Mb",
            "MethylMasteR–default bin size",
            "MethylMasteR–10kb",
            "MethylMasteR–100kb",
            "MethylMasteR–1Mb"
)


CCNE1 <- c(1:length(states))
CDKN2A <- c(1:length(states))
MYC <- c(1:length(states))
MYOCD <- c(1:length(states))
PTEN <- c(1:length(states))
RB1 <- c(1:length(states))
TP53 <- c(1:length(states))

output.df <- data.frame(states, CCNE1,
                        CDKN2A,
                        MYC,
                        MYOCD, 
                        PTEN ,
                        RB1 ,
                        TP53)

output.df.pearson <- output.df %>%
  separate(
    col = states,
    into = c("Tech", "Bin"),
    sep = "–"
  )

output.df.spearman <- output.df %>%
  separate(
    col = states,
    into = c("Tech", "Bin"),
    sep = "–"
  )


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


### ACCURACY ----
## LMS ----
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

### LM ----
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

## Alternative LM Accuracy Calculation
# Just use CNV Accuracy
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

## Alternative LMS Accuracy Calculation
# Just use CNV Accuracy
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

 ### Case Graphs ----
## LMS ----
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

## LM ----
for(c in LM_cases){
  for(b in bins){
    my_plot <- plot_cnv_segments(df = rbind(LMStt(STT = c, tech = "MethylMaster", bin = b), 
                                            LMStt(STT = c, tech = "Sesame", bin = b), 
                                            LMStt(STT = c, tech = "Conumee", bin = b)))
    dir.create(paste0("~/Work/Analysis/quarterCutoff/case_graphs/LM/",c), showWarnings = FALSE, recursive = TRUE)
    ggsave(filename = paste0("~/Work/Analysis/quarterCutoff/case_graphs/LM/",c,"/", c,"-",b,".pdf"), plot = my_plot, width = 24, height = 8)
  }
}

## NORMALS ----
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





### Correlations ----
# Building correlations for Accuracy to Genome Modified, Genomic Index to Accuracy, and Genomic Index to Genome Modified
## Accuracy - Genomic Index ----
gi <- read.csv(paste0(getwd(), "/quarterCutoff/Genomic_Index/genomicIndexLMS.csv")) %>% 
  dplyr::select(c(Case, Bin, type, gi)) %>%
  dplyr::filter(type != "SNP")
acc <- read.csv(paste0(getwd(), "/quarterCutoff/Accuracy/accuracyLMS.csv")) %>% 
  dplyr::select(c(Cases, Accuracy, Bin_Size, Technology)) %>% 
  dplyr::rename(Case = Cases) %>%
  dplyr::rename(Bin = Bin_Size) %>%
  dplyr::rename(type = Technology)
acc[acc == "Sesame"] <- "SeSAMe"
data <- dplyr::inner_join(gi, acc, by = c("Case", "Bin", "type"))

# Function that computes Spearman correlation on a resampled dataset
spearman_boot <- function(data, indices) {
  d <- data[indices, ]
  cor(d$Accuracy, d$gi, method = "spearman")
}

# Wrapper to run the bootstrap for one technology/bin subset
bootstrap_correlation <- function(df, n_boot = 2000) {
  set.seed(42)  # for reproducibility
  
  boot_result <- boot(data = df, statistic = spearman_boot, R = n_boot)
  
  ci <- boot.ci(boot_result, type = "bca")  
  tibble(
    correlation = boot_result$t0,
    ci_lower = ci$bca[4],
    ci_upper = ci$bca[5]
  )
}

correlations <- data %>%
  group_by(Bin, type) %>%
  group_modify(~ bootstrap_correlation(.x)) %>%
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
  main = "Correlation Heatmap: Accuracy vs. Genomic Index"
)

# No worthy data from correlation for LMs.
## Accuracy - Genomic Modified ----
genome_modified <- read.csv(paste0(getwd(), "/quarterCutoff/Genome_modified/genome_modifiedLMS.csv")) %>% 
  dplyr::select(c(Case, Bin, Tech, val)) %>%
  dplyr::filter(Tech != "SNP") %>%
  dplyr::mutate(Bin = as.numeric(Bin))
acc <- read.csv(paste0(getwd(), "/quarterCutoff/Accuracy/accuracyLMS.csv")) %>% 
  dplyr::select(c(Cases, Accuracy, Bin_Size, Technology)) %>% 
  dplyr::rename(Case = Cases) %>%
  dplyr::rename(Bin = Bin_Size) %>%
  dplyr::rename(Tech = Technology)
data <- dplyr::inner_join(genome_modified, acc, by = c("Case", "Bin", "Tech"))

# Function that computes Spearman correlation on a resampled dataset
spearman_boot <- function(data, indices) {
  d <- data[indices, ]
  cor(d$Accuracy, d$val, method = "spearman")
}

# Wrapper to run the bootstrap for one technology/bin subset
bootstrap_correlation <- function(df, n_boot = 2000) {
  set.seed(42)  # for reproducibility
  
  boot_result <- boot(data = df, statistic = spearman_boot, R = n_boot)
  
  ci <- boot.ci(boot_result, type = "bca")  
  tibble(
    correlation = boot_result$t0,
    ci_lower = ci$bca[4],
    ci_upper = ci$bca[5]
  )
}

correlations <- data %>%
  group_by(Bin, Tech) %>%
  group_modify(~ bootstrap_correlation(.x)) %>%
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
  main = "Correlation Heatmap: Accuracy vs. Genome Modified"
)
