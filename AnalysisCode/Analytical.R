library(GenomicRanges)
library(dplyr)
accuracyModel <- function(t1){
  hg19ranges <- c(chr1 = 249250621, 
                  chr2 = 243199373,
                  chr3 = 198022430,
                  chr4 = 191154276,
                  chr5 = 180915260,
                  chr6 = 171115067,
                  chr7 = 159138663,
                  chr8 = 146364022,
                  chr9 = 141213431,
                  chr10 = 135534747,
                  chr11 = 135006516,
                  chr12 = 133851895,
                  chr13 = 115169878,
                  chr14 = 107349540,
                  chr15 = 102531392,
                  chr16 = 90354753,
                  chr17 = 81195210,
                  chr18 = 78077248,
                  chr19 = 59128983,
                  chr20 = 63025520,
                  chr21 = 48129895,
                  chr22 = 51304566)
  
  pred_df <- t1 %>% filter(!(type == "SNP"))
  truth_df <- t1 %>% filter(type == "SNP")
  
  truth_gr <- GRanges(seqnames = truth_df$chrom,
                      ranges = IRanges(start = truth_df$loc.start, end = truth_df$loc.end),
                      CNV = truth_df$CNVStatus)
  
  pred_gr <- GRanges(seqnames = pred_df$chrom,
                     ranges = IRanges(start = pred_df$loc.start, end = pred_df$loc.end),
                     CNV = pred_df$CNVStatus)
  
  
  hits <- findOverlaps(pred_gr, truth_gr)
  # print(hits)
  
  overlap_ranges <- pintersect(pred_gr[queryHits(hits)], truth_gr[subjectHits(hits)])
  # print(overlap_ranges)
  
  overlap_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = mcols(pred_gr)$CNV[queryHits(hits)],
    truth_cnv = mcols(truth_gr)$CNV[subjectHits(hits)]
  )
  # print(overlap_df)
  
  weighted_matrix <- overlap_df %>%
    group_by(truth_cnv, pred_cnv) %>%
    summarise(total_bp = sum(width), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pred_cnv, values_from = total_bp, values_fill = 0)
  # print(weighted_matrix)
  
  TP_weighted <- sum(overlap_df$width[overlap_df$truth_cnv == overlap_df$pred_cnv])
  # print(TP_weighted)
  
  Total_weighted <- sum(overlap_df$width)
  # print(Total_weighted)
  
  accuracy_by_chr <- overlap_df %>%
    group_by(Chromosome) %>%
    summarise(
      TP_bp = sum(width[pred_cnv == truth_cnv]),
      Total_bp = sum(width),
      Accuracy = TP_bp / Total_bp
    )
  weighted_accuracy <- TP_weighted / Total_weighted
  list(accuracy_by_chr, weighted_accuracy)
}
fpCheck <- function(df) {
  library(GenomicRanges)
  library(dplyr)
  
  pred_df <- df %>% filter(!(type == "SNP"))
  pred_df$chrom <- as.character(pred_df$chrom)
  truth_df <- df %>% filter(type == "SNP")
  truth_df$chrom <- as.character(truth_df$chrom)
  
  # ---- Step 1: hg19 chromosome sizes (numeric chromosomes 1 to 22) ----
  hg19_chr_sizes <- data.frame(
    Chromosome = as.character(1:22),
    Genome_bp = c(249250621, 243199373, 198022430, 191154276, 180915260,
                  171115067, 159138663, 146364022, 141213431, 135534747,
                  135006516, 133851895, 115169878, 107349540, 102531392,
                  90354753, 81195210, 78077248, 59128983, 63025520,
                  48129895, 51304566)
  )
  
  # ---- Step 2: Convert to GRanges ----
  truth_gr <- GRanges(seqnames = truth_df$chrom,
                      ranges = IRanges(start = truth_df$loc.start, end = truth_df$loc.end),
                      CNV = truth_df$CNVStatus)
  
  pred_gr <- GRanges(seqnames = pred_df$chrom,
                     ranges = IRanges(start = pred_df$loc.start, end = pred_df$loc.end),
                     CNV = pred_df$CNVStatus)
  
  # ---- Step 3: Overlaps and TP ----
  hits <- findOverlaps(pred_gr, truth_gr)
  
  overlap_ranges <- pintersect(pred_gr[queryHits(hits)], truth_gr[subjectHits(hits)])
  pred_cnv <- mcols(pred_gr)$CNV[queryHits(hits)]
  truth_cnv <- mcols(truth_gr)$CNV[subjectHits(hits)]
  
  tp_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = pred_cnv,
    truth_cnv = truth_cnv
  ) %>%
    filter(pred_cnv == truth_cnv) %>%
    group_by(Chromosome) %>%
    summarise(TP_bp = sum(width), .groups = "drop")
  print(tp_df)
  
  # ---- Step 4: False Positives ----
  non_overlapping_pred <- pred_gr[-queryHits(hits)]
  fp_df <- as.data.frame(non_overlapping_pred) %>%
    filter(CNV != "Normal") %>%
    group_by(seqnames) %>%
    summarise(FP_bp = sum(width), .groups = "drop") %>% dplyr::rename(Chromosome = "seqnames")
  print(fp_df)
  
  # ---- Step 5: False Negatives ----
  non_overlapping_truth <- truth_gr[-subjectHits(hits)]
  fn_df <- as.data.frame(non_overlapping_truth) %>%
    filter(CNV != "Normal") %>%
    group_by(seqnames) %>%
    summarise(FN_bp = sum(width), .groups = "drop") %>% dplyr::rename(Chromosome = "seqnames")
  print(fn_df)
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
    arrange(as.numeric(Chromosome))
  
  accuracy <- sum(final_eval$Accuracy) / 22
  
  return(list(final_eval, accuracy))
}
pccB <- function(dset){
  truth_df <- dset %>% filter(type == "SNP")
  test_df <- dset %>% filter(!(type == "SNP"))
  
  truth_gr <- GRanges(seqnames = truth_df$chrom,
                      ranges = IRanges(start = truth_df$loc.start, end = truth_df$loc.end),
                      CNV = truth_df$seg.mean)
  
  pred_gr <- GRanges(seqnames = test_df$chrom,
                     ranges = IRanges(start = test_df$loc.start, end = test_df$loc.end),
                     CNV = test_df$seg.mean)
  
  
  hits <- findOverlaps(pred_gr, truth_gr)
  
  overlap_ranges <- pintersect(pred_gr[queryHits(hits)], truth_gr[subjectHits(hits)])
  
  overlap_widths <- width(overlap_ranges)
  
  overlap_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = overlap_ranges,
    pred_cnv = mcols(pred_gr)$CNV[queryHits(hits)],
    truth_cnv = mcols(truth_gr)$CNV[subjectHits(hits)]
  ) %>% dplyr::select(-c(width.seqnames, width.strand, width.CNV, width.hit))
  
  # overlap_df 
  cor(overlap_df$pred_cnv, overlap_df$truth_cnv, method = "pearson")
  
}

LMMethylCalling <- function(df){
  
} # FOCAL CNV CALLING METHYLATION DATA
LMSMethylCalling <- function(df){} # WIDE CNV CALLING METHYLATION DATA
LMcallingSNP <- function(Sample) {} # Focal CNV SNP
LMScallingSNP <- function(Sample) {
  commonCSVFolder <- "~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV"
  readEvents <- read.csv(paste0(commonCSVFolder, "/", Sample, "_events.csv")) %>%
    dplyr::select("Chromosome.Region", "Probe.Median", "Event", "Length") %>%
    separate_wider_delim(
      cols = "Chromosome.Region",
      delim = ":",
      names = c("chrom", "locs")
    ) %>%
    separate_wider_delim(
      cols = "locs",
      delim = "-",
      names = c("loc.start", "loc.end")
    ) %>% dplyr::mutate(loc.start = as.numeric(str_remove_all(loc.start, ",")), 
                        loc.end = as.numeric(str_remove_all(loc.end, ",")),
                        chrom = as.numeric(str_remove_all(chrom, "chr"))) %>%
    dplyr::mutate(CNVStatus = ifelse(Event == "CN Loss", "Deletion", "Amplification")) %>%
    dplyr::select(-c("Event")) %>% dplyr::rename(seg.mean = Probe.Median) %>% 
    dplyr::filter(!is.na(chrom))
  
  # Filtering for wider SNPs (200 probes as provided in planning document)
  # Developed via oncoscanR for sourcing the OncoScan manifest. - Write w/ internet
  # OncoScan.na33.r4 is also of interest for this bit
  
  
  
} # Wide CNV SNP


# Genomic Index = A^2/22 (A being gains & losses sum.)
GenomicIndex <- function(dfCH3, sw){
  if(sw == F){
    dfSNP <- dfCH3 %>% dplyr::filter(type == "SNP") %>%
      dplyr::filter(CNVStatus != "Normal")
    CNVCount <- nrow(dfSNP)
    ret <- (CNVCount^2)/22
    return(ret)
  } else {
  
  # Filter out for clinically significant segments %>% filter(CNVStatus != Normal)
  dfCH3f <- dfCH3 %>% dplyr::filter(type == "Conumee" || type == "MethylMasteR" || type == "SeSAMe") %>%
    dplyr::filter(CNVStatus != "Normal")
  print(nrow(dfCH3))
  # Count
  CNVCount <- nrow(dfCH3f)
  
  
  # Return value
  ret <- (CNVCount^2)/22
  print(CNVCount)
  return(ret) 
  }
}

# Genomic Index runner.
correlative <- read_csv("LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
Def <- c()
kb10<- c()
kb100 <- c()
mb1 <- c()
trt <- c()
for(cr in correlative){
  Def <- c(Def, GenomicIndex(labLMSProc(cr, "Conumee", 50000), sw = T))
  kb10 <- c(kb10, GenomicIndex(labLMSProc(cr, "Conumee", 10000), sw = T))
  kb100 <- c(kb100, GenomicIndex(labLMSProc(cr, "Conumee", 1e+05), sw = T))
  mb1 <- c(mb1, GenomicIndex(labLMSProc(cr, "Conumee", 1e+06), sw = T))
  trt <- c(trt, GenomicIndex(labLMSProc(cr, "Conumee", 50000), sw = F))
}
GIndex <- data.frame(STT = correlative, ConumeeDefault = Def, Conumee10kb = kb10, Conumee100kb = kb100, Conumee1mb = mb1, BaseTruth = trt)

Def <- c()
kb10<- c()
kb100 <- c()
mb1 <- c()
for(cr in correlative){
  Def <- c(Def, GenomicIndex(labLMSProc(cr, "Sesame", 50000), sw = T))
  kb10 <- c(kb10, GenomicIndex(labLMSProc(cr, "Sesame", 10000), sw = T))
  kb100 <- c(kb100, GenomicIndex(labLMSProc(cr, "Sesame", 1e+05), sw = T))
  mb1 <- c(mb1, GenomicIndex(labLMSProc(cr, "Sesame", 1e+06), sw = T))
}
GIndex <- GIndex %>% dplyr::mutate(SesameDefault = Def, Sesame10kb = kb10, 
                                   Sesame100kb = kb100, Sesame1mb = mb1)
