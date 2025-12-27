# Rough diploid recentering algorithm that applies averaging across all detected 
# segments.

# genomicLoc is represented as a dataframe of one row with labels which 
# give the chromosome, the start and the end.
require(GenomicRanges)
require(tidyverse)
genomicLoc <- data.frame(chromosome = 1, start = 160000000, end = 200000000)
dipRecenterAndCorrect <- function(cbs, genomicLoc){
  cbsF <- read_delim(cbs)
  std <- cbsF %>% dplyr::filter(Chr == genomicLoc$chromosome[1])
  std_gr <- makeGRangesFromDataFrame(std, keep.extra.columns = TRUE)
  gLoc_gr <- makeGRangesFromDataFrame(genomicLoc)
  ovl <- findOverlaps(std_gr, gLoc_gr)
  overlaps_intervals <- pintersect(std_gr[queryHits(ovl)], gLoc_gr[subjectHits(ovl)])
  finalDf <- as.data.frame(overlaps_intervals) %>% select(c("seqnames", "start", "end", "Log2Ratio")) %>%
    dplyr::rename(chr = "seqnames", seg.mean = "Log2Ratio", loc.start = "start", loc.end = "end")
  stdMed <- median(finalDf$seg.mean)
  cbsF <- cbsF %>% dplyr::mutate(Log2Ratio = Log2Ratio - stdMed) %>% dplyr::filter(Chr != 23 && Chr != 24)
  # write_tsv(cbsF, file = cbs)
  cbsF
}

# Sample tester to compare and contrast the two outputs.
testCorr <- function(sampleName){
  tst <- dipRecenterAndCorrect(paste0("~/Documents/Analysis/LabData/ReprocessedLMS/",sampleName,"/ASCAT/L2R/",sampleName,".NoCut.cbs"),
                               genomicLoc = genomicLoc)
  tstSEQ <- dipRecenterAndCorrect(paste0("~/Documents/Analysis/LabData/ReprocessedLMS/",sampleName,"/SEQUENZA/L2R/",sampleName,".NoCut.cbs"),
                                  genomicLoc = genomicLoc)
  FASST2 <- read_tsv(paste0("~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/",sampleName,"_events.txt"), 
                     skip = 36) 
  FASST2 <- rbind(FASST2 %>% dplyr::filter(Event == "CN Loss"), 
                  FASST2 %>% dplyr::filter(Event == "CN Gain")) %>%
    dplyr::select(c("Probe Median", "Event", "Chromosome Region"))
  FASST2 <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(FASST2, `Chromosome Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
                                        names = c("loc.start", "loc.end"), delim = "-")
  FASST2$loc.start <- as.numeric(str_replace_all(FASST2$loc.start, ",", ""))
  FASST2$loc.end <- as.numeric(str_replace_all(FASST2$loc.end, ",", ""))
  FASST2$chrom <- as.numeric(str_replace_all(FASST2$chrom, "chr", ""))
  FASST2 <- drop_na(FASST2) 
  tst_gr <- makeGRangesFromDataFrame(tst, keep.extra.columns = TRUE)
  FASST2_gr <- makeGRangesFromDataFrame(FASST2, keep.extra.columns = TRUE)
  ovl <- findOverlaps(tst_gr, FASST2_gr)
  gr1_overlap <- tst_gr[queryHits(ovl)]
  gr2_overlap <- FASST2_gr[subjectHits(ovl)]
  mcols(gr1_overlap) <- cbind(mcols(gr1_overlap), mcols(gr2_overlap))
  comp_df <- as.data.frame(gr1_overlap) %>% 
    select(c("seqnames", "start", "end", "Log2Ratio", "Probe.Median")) %>% 
    dplyr::rename(chrom = "seqnames", loc.start = "start", loc.end = "end", 
                  ASCAT = "Log2Ratio", FASST2 = "Probe.Median") %>%
    dplyr::mutate(FASST2 = as.numeric(FASST2))
  cor(comp_df$ASCAT, comp_df$FASST2, method = "pearson")
}
testCorrSEQUENZA <- function(sampleName){
  tst <- dipRecenterAndCorrect(paste0("~/Documents/Analysis/LabData/ReprocessedLMS/",sampleName,"/SEQUENZA/L2R/",sampleName,".NoCut.cbs"),
                                  genomicLoc = genomicLoc)
  FASST2 <- read_tsv(paste0("~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/",sampleName,"_events.txt"), 
                     skip = 36) 
  FASST2 <- rbind(FASST2 %>% dplyr::filter(Event == "CN Loss"), 
                  FASST2 %>% dplyr::filter(Event == "CN Gain")) %>%
    dplyr::select(c("Probe Median", "Event", "Chromosome Region"))
  FASST2 <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(FASST2, `Chromosome Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
                                        names = c("loc.start", "loc.end"), delim = "-")
  FASST2$loc.start <- as.numeric(str_replace_all(FASST2$loc.start, ",", ""))
  FASST2$loc.end <- as.numeric(str_replace_all(FASST2$loc.end, ",", ""))
  FASST2$chrom <- as.numeric(str_replace_all(FASST2$chrom, "chr", ""))
  FASST2 <- drop_na(FASST2) 
  tst_gr <- makeGRangesFromDataFrame(tst, keep.extra.columns = TRUE)
  FASST2_gr <- makeGRangesFromDataFrame(FASST2, keep.extra.columns = TRUE)
  ovl <- findOverlaps(tst_gr, FASST2_gr)
  gr1_overlap <- tst_gr[queryHits(ovl)]
  gr2_overlap <- FASST2_gr[subjectHits(ovl)]
  mcols(gr1_overlap) <- cbind(mcols(gr1_overlap), mcols(gr2_overlap))
  comp_df <- as.data.frame(gr1_overlap) %>% 
    select(c("seqnames", "start", "end", "Log2Ratio", "Probe.Median")) %>% 
    dplyr::rename(chrom = "seqnames", loc.start = "start", loc.end = "end", 
                  ASCAT = "Log2Ratio", FASST2 = "Probe.Median") %>%
    dplyr::mutate(FASST2 = as.numeric(FASST2))
  cor(comp_df$ASCAT, comp_df$FASST2, method = "pearson")
  
}
spExtra <- read_excel("~/Documents/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/design_13LMS_CNVs_other_info_12Aug2025.xlsx")
correlations <- c()
correlationsSEQ <- c()
for(sp in spExtra$Sample){
  correlations <- c(testCorr(sp), correlations)
  correlationsSEQ <- c(testCorrSEQUENZA(sp), correlationsSEQ)
}
spExtra <- spExtra %>% dplyr::mutate(correlation_to_ASCAT = correlations, correlation_to_SEQUENZA = correlationsSEQ) %>%
  dplyr::select("Sample", "STT", "Region used for              diploid re-centering", "correlation_to_ASCAT", "correlation_to_SEQUENZA")
write.csv(spExtra, file = "~/Documents/Analysis/LabData/ReprocessedLMS/correlation_w_FASST2.csv")
