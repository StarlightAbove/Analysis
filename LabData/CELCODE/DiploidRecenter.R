# Rough diploid recentering algorithm that applies averaging across all detected 
# segments.

# genomicLoc is represented as a dataframe of one row with labels which 
# give the chromosome, the start and the end.
require(GenomicRanges)
require(tidyverse)
library(stringr)
genomicLoc <- data.frame(chromosome = 1, start = 160000000, end = 200000000)
dipRecenterAndCorrect <- function(rds, genomicLoc = NULL){
  rdsR <- readRDS(rds)
  if(is.null(genomicLoc)){
    # Plan:
    # .RDS value first passed, raw.
    genome.pkg <- rdsR$meta$basic$genome.pkg
    BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
    cs <- chromobjector(BSg.obj)
    l2r.chr <- unname(unlist(cs$chrom2chr[as.character(rdsR$data$SNPpos$chrs)]))
    baf.value <- data.frame(Chr = l2r.chr,
    Start = as.integer(rdsR$data$SNPpos$pos),
    End = as.integer(rdsR$data$SNPpos$pos),
    Value = rdsR$data$Tumor_BAF[,1],
    stringsAsFactors = FALSE)
    print(baf.value)
    #
    # All of this hijacked from EaCoN.
    #
    # Now, read the B-allele frequencies and via their locations, look for homozygous (0.0) as candidate regions.
    # Determine erroneous significance by matching that to anomalous log2 ratio by passing SEG.SEGMENTER.RDS and filter.
    # Check against Gene database to determine whether there are any genes in the candidate region and filter.
    # Print candidates and ask for input or select best.
    # Subtract median across all log-2 ratios in the SEG RDS.
    
  }
  #std <- cbsF %>% dplyr::filter(Chr == genomicLoc$chromosome[1])
  #std_gr <- makeGRangesFromDataFrame(std, keep.extra.columns = TRUE)
  #gLoc_gr <- makeGRangesFromDataFrame(genomicLoc)
  #ovl <- findOverlaps(std_gr, gLoc_gr)
  #overlaps_intervals <- pintersect(std_gr[queryHits(ovl)], gLoc_gr[subjectHits(ovl)])
  #finalDf <- as.data.frame(overlaps_intervals) %>% select(c("seqnames", "start", "end", "Log2Ratio")) %>%
  #  dplyr::rename(chr = "seqnames", seg.mean = "Log2Ratio", loc.start = "start", loc.end = "end")
  #stdMed <- median(finalDf$seg.mean)
  #cbsF <- cbsF %>% dplyr::mutate(Log2Ratio = Log2Ratio - stdMed) %>% dplyr::filter(Chr != 23 && Chr != 24)
  # write_tsv(cbsF, file = cbs)
  #cbsF
}

dipRecenterAndCorrect("~/Work/Analysis/LabData/ReprocessedLMS/10_CAD_9355_LMS_Rec3/10_CAD_9355_LMS_Rec3_OncoScan_hg19_processed.RDS")
