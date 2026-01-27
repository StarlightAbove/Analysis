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
    
    my_tibble <- data.frame(Bafseg = rdsR$data$Tumor_BAF_segmented[[1]])
    keeps <- rownames(my_tibble)
    correl <- rdsR$data$SNPpos %>% dplyr::filter(rownames(rdsR$data$SNPpos) %in% keeps) # Contains location data for each SNP.
    
    
    baf.value <- data.frame(Chr = l2r.chr,
    Start = as.integer(rdsR$data$SNPpos$pos),
    End = as.integer(rdsR$data$SNPpos$pos),
    Value = rdsR$data$Tumor_BAF[,1],
    stringsAsFactors = FALSE)
    Segmented_BAF_loc <- merge(correl, my_tibble, by = "row.names")
    Segmented_BAF_loc$chrs <- as.numeric(str_remove_all(Segmented_BAF_loc$chrs, "chr"))
    Segmented_BAF_loc <- Segmented_BAF_loc %>% dplyr::filter(!is.na(chrs)) %>% 
      dplyr::arrange(chrs, pos) %>% dplyr::mutate(group_id = consecutive_id(chrs, Bafseg)) 
    SegmentLoc <- Segmented_BAF_loc %>% group_by(group_id) %>%
      summarise(
        chrom = mean(chrs),
        Start = min(pos, na.rm = TRUE),
        End = max(pos, na.rm = TRUE),
        BAF = median(Bafseg),
        .groups = 'drop'
      ) %>% dplyr::filter(BAF == 0 | BAF == 0.5 | BAF == 1.0) %>% dplyr::select(-c("group_id"))
    # We've gotten the candidate regions based on BAF.
    
    l2s <- rdsR[["cbs"]][["cut"]] %>% dplyr::select(c("Chr", "Start", "End", "Log2Ratio"))
    l2s_gr <- GenomicRanges::makeGRangesFromDataFrame(l2s, keep.extra.columns = T)
    segLoc_gr <- GenomicRanges::makeGRangesFromDataFrame(SegmentLoc, keep.extra.columns = T)
    sub <- findOverlaps(l2s_gr, segLoc_gr)
    overlaps_intervals <- pintersect(l2s_gr[queryHits(sub)], segLoc_gr[subjectHits(sub)])
    overlaps_intervals$BAF <- 0.5
    overlaps_intervals <- as.data.frame(overlaps_intervals) %>% dplyr::select(-c("hit")) %>%
      dplyr::filter(!(Log2Ratio < 0.2 & Log2Ratio > -0.2) ) %>% dplyr::arrange(desc(width), desc(Log2Ratio))
    
    # Question to ask, which corrective factor to use?
    options <- c(paste0("chr",overlaps_intervals$seqnames,": ",overlaps_intervals$start, "-", overlaps_intervals$end, " (", overlaps_intervals$Log2Ratio,")"))
    selected <- menu(options, title = "Select a correction location: ")
    corrector <- overlaps_intervals[selected,]$Log2Ratio[1]
    rdsR$cbs$cut$Log2Ratio <- rdsR$cbs$cut$Log2Ratio - corrector
    write_rds(x = rdsR, file = paste0("~/Work/Analysis/LabData/ReprocessedLMS/DiploidCorr/",colnames(rdsR$data$Tumor_BAF),".ASCAT.RDS"))
    
    #print(overlaps_intervals[selected])
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

abc <- dipRecenterAndCorrect("~/Work/Analysis/LabData/ReprocessedLMS/8_KAN_9337_LMS_Rec/ASCAT/L2R/8_KAN_9337_LMS_Rec.SEG.ASCAT.RDS")
