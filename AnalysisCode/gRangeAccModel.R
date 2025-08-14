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
  
  overlap_ranges <- pintersect(pred_gr[queryHits(hits)], truth_gr[subjectHits(hits)])
  
  overlap_widths <- width(overlap_ranges)
  
  overlap_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = mcols(pred_gr)$CNV[queryHits(hits)],
    truth_cnv = mcols(truth_gr)$CNV[subjectHits(hits)]
  )
  
  weighted_matrix <- overlap_df %>%
    group_by(truth_cnv, pred_cnv) %>%
    summarise(total_bp = sum(width), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pred_cnv, values_from = total_bp, values_fill = 0)
  
  accuracy_by_chr <- overlap_df %>%
    group_by(Chromosome) %>%
    summarise(
      TP_bp = sum(width[pred_cnv == truth_cnv]),
      Total_bp = sum(width),
      Accuracy = TP_bp / Total_bp
    )
  
  TP_weighted <- sum(overlap_df$width[overlap_df$truth_cnv == overlap_df$pred_cnv])
  
  Total_weighted <- sum(overlap_df$width)
  
  weighted_accuracy <- TP_weighted / Total_weighted
  list(accuracy_by_chr, weighted_accuracy)
}

# Issue log:
# - LM Intersections aren't as accurate as it does not account for correct normal predictions.
# Solution: 
# - Presume that normal-normal matching in empty chromosomes means the chromosome has an accuracy of 1.