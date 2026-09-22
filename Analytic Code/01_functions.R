### ============================================================
### 1. FUNCTIONS
### ============================================================
# All shared helper functions. Pure definitions, no side effects.

# All helper functions used by the analysis sections below. Grouped
# roughly as: (a) per-case data loaders/mergers for each cohort,
# (b) genomic-index / genome-modified / accuracy metric functions,
# (c) gene annotation + plotting helpers.

#### 1a. Per-case data loaders (read raw caller output, standardize columns) ----
# Load & standardize one LMS case's CNV segments for a given caller
# ("MethylMaster", "Conumee", or "Sesame"), combined with its matched
# SNP-array (ChAS) segments as a common ground-truth reference.
labLMSProc <- function(STTq, Technology, binSize){
  
  
  # Correlate by STT information between methylation Sentrix and SNP data.
  correlationSheet <- read.csv(paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/correlative.csv")) %>% filter(STT == STTq)
  
  cnvMatch <- read_delim(paste0(getwd(), "/LabData/LMS_SNP_EPIC_array_data/ChAS/ChAS_data_01Feb2026/ChAS_LMS_Probe_and_segment_level_data_01Feb2026/STT", 
                                STTq, "_Recentered_Segment_level_data_01Feb2026.segment.txt"))
  
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
# LM-cohort equivalent of labLMSProc(): loads the SNP-array segments for
# one LM case plus that case's segments from a single named caller ("tech").
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
#### 1b. Normal-control loader ----
# Load & standardize one normal-tissue control's CNV segments for a
# given caller (there is no SNP ground truth for normals).
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


# Summarize each caller's median segment log2-ratio across a set of
# normal controls (IDs) at a given bin size, as a sanity check that
# normals cluster near zero / within an acceptable cutoff.
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
#### 1c. Metric functions (genomic index, genome modified, accuracy) ----

# Genomic Index = (number of non-normal CNV segments)^2 / (number of
# distinct chromosomes with a call). Computed either directly from the
# raw SNP-array file (when sw == FALSE) or from an already-combined
# data frame (dfCH3) for a methylation caller.
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


# Build the row-level (non-aggregated) table of non-normal CNV segments
# for one case/bin/technology, combining the caller's calls with the
# matching SNP-array calls. Used as the input to the group_by/summarize
# step that actually computes the Genomic Index (see section 4).
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


# % of the (22-autosome) hg19 genome covered by non-normal CNV segments.
GenomeModified <- function(dfCH3){
  df <- dfCH3 %>% dplyr::filter(type != "SNP", CNVStatus != "Normal") %>%
    dplyr::mutate(width = loc.end - loc.start) 
  modifiedWidth <- sum(df$width)
  hg19_info <- getChromInfoFromUCSC("hg19") %>% dplyr::filter(assembled == "TRUE")
  hg19_info <- hg19_info[1:22,]
  hg19_total <- sum(hg19_info$size)
  
  return((modifiedWidth/hg19_total) * 100)
  
}


# Per-chromosome true/false positive/negative base-pair accounting for a
# methylation caller against the SNP-array ground truth, plus an overall
# per-case accuracy summary. Returns list(per_chromosome_table, median_accuracy).
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
#### 1d. Gene annotation & plotting helpers ----

# Look up genomic coordinates for a set of gene symbols (via org.Hs.eg.db
# + the hg19 TxDb) and, if a combined CNV data frame (db) is supplied,
# overlap those coordinates onto it to pull out (or interpolate, for
# genes that fall between segments) the log2 ratio at each gene per
# caller. Data-only version of geneGen() below (no heatmap produced).
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


# Whole-genome CNV segment plot: log2 ratio vs. cumulative genomic
# position, one colored line per caller, chromosome boundaries marked,
# with gene labels overlaid where the input df contains "Gene_*" rows.
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

# Like NoGraphGeneGen(), but also builds a per-gene x per-caller heatmap
# of log2 ratios. Returns list(pheatmap_object, overlapping_regions_df).
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

# Combine geneGen()'s heatmap with plot_cnv_segments()'s genome-wide
# track into a single side-by-side (patchwork) figure for one case.
geneAnno <- function(Gene, db = NULL){
  reference <- geneGen(Gene = Gene, db = db)
  db <- rbind(db, reference[[2]])
  heatmap <- reference[[1]]
  pheatmap_ggplot <- as.ggplot(heatmap$gtable)
  plot <- plot_cnv_segments(db)
  return(plot + pheatmap_ggplot + plot_layout(widths = unit(c(20,8), c("null","null"))))
}
