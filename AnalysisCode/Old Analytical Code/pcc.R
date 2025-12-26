# Pearson correlation coefficient processor, by block and by whole dataset.
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

# Lab LMS correlation.
correlative <- read_csv("LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
casesList <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", 
               "Sample 6", "Sample 7", "Sample 8", "Sample 9", "Sample 10", 
               "Sample 11", "Sample 12", "Sample 13")
ses <- c(pccB(labLMSProc(correlative[1], "Sesame")), 
         pccB(labLMSProc(correlative[2], "Sesame")), 
         pccB(labLMSProc(correlative[3], "Sesame")), 
         pccB(labLMSProc(correlative[4], "Sesame")),
         pccB(labLMSProc(correlative[5], "Sesame")),
         pccB(labLMSProc(correlative[6], "Sesame")),
         pccB(labLMSProc(correlative[7], "Sesame")),
         pccB(labLMSProc(correlative[8], "Sesame")),
         pccB(labLMSProc(correlative[9], "Sesame")),
         pccB(labLMSProc(correlative[10], "Sesame")),
         pccB(labLMSProc(correlative[11], "Sesame")),
         pccB(labLMSProc(correlative[12], "Sesame")),
         pccB(labLMSProc(correlative[13], "Sesame")))

con <- c(pccB(labLMSProc(correlative[1], "Conumee")), 
         pccB(labLMSProc(correlative[2], "Conumee")), 
         pccB(labLMSProc(correlative[3], "Conumee")), 
         pccB(labLMSProc(correlative[4], "Conumee")),
         pccB(labLMSProc(correlative[5], "Conumee")),
         pccB(labLMSProc(correlative[6], "Conumee")),
         pccB(labLMSProc(correlative[7], "Conumee")),
         pccB(labLMSProc(correlative[8], "Conumee")),
         pccB(labLMSProc(correlative[9], "Conumee")),
         pccB(labLMSProc(correlative[10], "Conumee")),
         pccB(labLMSProc(correlative[11], "Conumee")),
         pccB(labLMSProc(correlative[12], "Conumee")),
         pccB(labLMSProc(correlative[13], "Conumee")))

methyl <- c(pccB(labLMSProc(correlative[1], "MethylMaster")), 
            pccB(labLMSProc(correlative[2], "MethylMaster")), 
            pccB(labLMSProc(correlative[3], "MethylMaster")), 
            pccB(labLMSProc(correlative[4], "MethylMaster")),
            pccB(labLMSProc(correlative[5], "MethylMaster")),
            pccB(labLMSProc(correlative[6], "MethylMaster")),
            pccB(labLMSProc(correlative[7], "MethylMaster")),
            pccB(labLMSProc(correlative[8], "MethylMaster")),
            pccB(labLMSProc(correlative[9], "MethylMaster")),
            pccB(labLMSProc(correlative[10], "MethylMaster")),
            pccB(labLMSProc(correlative[11], "MethylMaster")),
            pccB(labLMSProc(correlative[12], "MethylMaster")),
            pccB(labLMSProc(correlative[13], "MethylMaster")))
cMatrix <- data.frame(Cases = casesList, STT = correlative, 
                      Conumee = con, SeSAMe = ses, MethylMaster = methyl)
pltM <- pivot_longer(cMatrix, names_to = "Technology", values_to = "Correlation", cols = c("Conumee", "SeSAMe", "MethylMaster"))
ggplot(pltM, aes(fill=Technology, y=Correlation, x=Cases)) + 
  geom_bar(position="dodge", stat="identity")

# LMS Correlation.
barcodes <- as.list(read.csv("./cases/SingleFileCase/SampleSheet.csv") %>% select("caseNames"))[[1]]
casesList <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", 
               "Sample 6", "Sample 7", "Sample 8", "Sample 9", "Sample 10")
ses <- c(pccB(softOutputProcessingLMS(barcodes[1])), 
         pccB(softOutputProcessingLMS(barcodes[2])), 
         pccB(softOutputProcessingLMS(barcodes[3])), 
         pccB(softOutputProcessingLMS(barcodes[4])),
         pccB(softOutputProcessingLMS(barcodes[5])),
         pccB(softOutputProcessingLMS(barcodes[6])),
         pccB(softOutputProcessingLMS(barcodes[7])),
         pccB(softOutputProcessingLMS(barcodes[8])),
         pccB(softOutputProcessingLMS(barcodes[9])),
         pccB(softOutputProcessingLMS(barcodes[10])))

con <- c(pccB(matchingAlgo(barcodes[1], "Conumee", "LMS")), 
         pccB(matchingAlgo(barcodes[2], "Conumee", "LMS")), 
         pccB(matchingAlgo(barcodes[3], "Conumee", "LMS")), 
         pccB(matchingAlgo(barcodes[4], "Conumee", "LMS")),
         pccB(matchingAlgo(barcodes[5], "Conumee", "LMS")),
         pccB(matchingAlgo(barcodes[6], "Conumee", "LMS")),
         pccB(matchingAlgo(barcodes[7], "Conumee", "LMS")),
         pccB(matchingAlgo(barcodes[8], "Conumee", "LMS")),
         pccB(matchingAlgo(barcodes[9], "Conumee", "LMS")),
         pccB(matchingAlgo(barcodes[10], "Conumee", "LMS")))

methyl <- c(pccB(lmsProcessing(barcodes[1])), 
            pccB(lmsProcessing(barcodes[2])), 
            pccB(lmsProcessing(barcodes[3])), 
            pccB(lmsProcessing(barcodes[4])),
            pccB(lmsProcessing(barcodes[5])),
            pccB(lmsProcessing(barcodes[6])),
            pccB(lmsProcessing(barcodes[7])),
            pccB(lmsProcessing(barcodes[8])),
            pccB(lmsProcessing(barcodes[9])),
            pccB(lmsProcessing(barcodes[10])))
            
cMatrix <- data.frame(Cases = casesList, STT = barcodes, 
                      Conumee = con, SeSAMe = ses, MethylMaster = methyl)
pltM <- pivot_longer(cMatrix, names_to = "Technology", values_to = "Correlation", cols = c("Conumee", "SeSAMe", "MethylMaster"))
ggplot(pltM, aes(fill=Technology, y=Correlation, x=Cases)) + 
  geom_bar(position="dodge", stat="identity")

# LM Correlation.
barcodesCM <- paste0("./Outputs/Conumee/LMData/", list.files("./Outputs/Conumee/LMData/", pattern = "\\.csv$") ) 
barcodesSMLM <- list.files("./Outputs/SeSAMe/LM/") 
barcodesSMLM <- paste0("./Outputs/SeSAMe/LM/", barcodesSMLM[!(grepl(pattern = "bins_", barcodesSMLM))])
barcodesMMLM <- paste0("./Outputs/MethylMaster/LM/", list.files("./Outputs/MethylMaster/LM", pattern = "autocorrected_regions.csv", recursive = T))
casesList <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", 
               "Sample 6", "Sample 7", "Sample 8", "Sample 9")
ses <- c(pccB(softOutputProcessingLM(barcodesSMLM[1])), 
         pccB(softOutputProcessingLM(barcodesSMLM[2])), 
         pccB(softOutputProcessingLM(barcodesSMLM[3])), 
         pccB(softOutputProcessingLM(barcodesSMLM[4])),
         pccB(softOutputProcessingLM(barcodesSMLM[5])),
         pccB(softOutputProcessingLM(barcodesSMLM[6])),
         pccB(softOutputProcessingLM(barcodesSMLM[7])),
         pccB(softOutputProcessingLM(barcodesSMLM[8])),
         pccB(softOutputProcessingLM(barcodesSMLM[9])))

con <- c(pccB(matchingAlgoLM(barcodesCM[1])), 
         pccB(matchingAlgoLM(barcodesCM[1])), 
         pccB(matchingAlgoLM(barcodesCM[1])), 
         pccB(matchingAlgoLM(barcodesCM[1])),
         pccB(matchingAlgoLM(barcodesCM[1])),
         pccB(matchingAlgoLM(barcodesCM[1])),
         pccB(matchingAlgoLM(barcodesCM[1])),
         pccB(matchingAlgoLM(barcodesCM[1])),
         pccB(matchingAlgoLM(barcodesCM[1])))

methyl <- c(pccB(lmProcessing(barcodesMMLM[1])), 
            pccB(lmProcessing(barcodesMMLM[2])), 
            pccB(lmProcessing(barcodesMMLM[3])), 
            pccB(lmProcessing(barcodesMMLM[4])),
            pccB(lmProcessing(barcodesMMLM[5])),
            pccB(lmProcessing(barcodesMMLM[6])),
            pccB(lmProcessing(barcodesMMLM[7])),
            pccB(lmProcessing(barcodesMMLM[8])),
            pccB(lmProcessing(barcodesMMLM[9])))

cMatrix <- data.frame(Cases = casesList, STT = barcodesCM, 
                      Conumee = con, SeSAMe = ses, MethylMaster = methyl)
pltM <- pivot_longer(cMatrix, names_to = "Technology", values_to = "Correlation", cols = c("Conumee", "SeSAMe", "MethylMaster"))
ggplot(pltM, aes(fill=Technology, y=Correlation, x=Cases)) + 
  geom_bar(position="dodge", stat="identity")


