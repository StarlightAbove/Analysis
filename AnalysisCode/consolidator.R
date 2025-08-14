# Consolidates data from cases into a single dataframe for easy analysis.
# Generalized barcodes for LMS to correlate case names.
barcodes <- as.list(read.csv("./cases/SingleFileCase/SampleSheet.csv") %>% select("caseNames"))[[1]]
# Conumee. - LMS
case1 <- accuracyModel(matchingAlgo(barcodes[1], "Conumee", "LMS"))
case2 <- accuracyModel(matchingAlgo(barcodes[2], "Conumee", "LMS"))
case3 <- accuracyModel(matchingAlgo(barcodes[3], "Conumee", "LMS"))
case4 <- accuracyModel(matchingAlgo(barcodes[4], "Conumee", "LMS"))
case5 <- accuracyModel(matchingAlgo(barcodes[5], "Conumee", "LMS"))
case6 <- accuracyModel(matchingAlgo(barcodes[6], "Conumee", "LMS"))
case7 <- accuracyModel(matchingAlgo(barcodes[7], "Conumee", "LMS"))
case8 <- accuracyModel(matchingAlgo(barcodes[8], "Conumee", "LMS"))
case9 <- accuracyModel(matchingAlgo(barcodes[9], "Conumee", "LMS"))
case10 <- accuracyModel(matchingAlgo(barcodes[10], "Conumee", "LMS"))

# SeSAMe - LMS
seSAMe1 <- accuracyModel(softOutputProcessingLMS(barcodes[1]))
seSAMe2 <- accuracyModel(softOutputProcessingLMS(barcodes[2]))
seSAMe3 <- accuracyModel(softOutputProcessingLMS(barcodes[3]))
seSAMe4 <- accuracyModel(softOutputProcessingLMS(barcodes[4]))
seSAMe5 <- accuracyModel(softOutputProcessingLMS(barcodes[5])) 
seSAMe6 <- accuracyModel(softOutputProcessingLMS(barcodes[6])) 
seSAMe7 <- accuracyModel(softOutputProcessingLMS(barcodes[7]))
seSAMe8 <- accuracyModel(softOutputProcessingLMS(barcodes[8]))
seSAMe9 <- accuracyModel(softOutputProcessingLMS(barcodes[9]))
seSAMe10 <- accuracyModel(softOutputProcessingLMS(barcodes[10]))

# MethylMasteR - LMS
mmstr1 <- accuracyModel(lmsProcessing(barcodes[1]))
mmstr2 <- accuracyModel(lmsProcessing(barcodes[2]))
mmstr3 <- accuracyModel(lmsProcessing(barcodes[3]))
mmstr4 <- accuracyModel(lmsProcessing(barcodes[4]))
mmstr5 <- accuracyModel(lmsProcessing(barcodes[5]))
mmstr6 <- accuracyModel(lmsProcessing(barcodes[6]))
mmstr7 <- accuracyModel(lmsProcessing(barcodes[7])) 
mmstr8 <- accuracyModel(lmsProcessing(barcodes[8])) 
mmstr9 <- accuracyModel(lmsProcessing(barcodes[9]))
mmstr10 <- accuracyModel(lmsProcessing(barcodes[10]))


# LMS Consolidation.
# Consolidation Design: Key - Cases, Chromosomes; Columns - Software
ses <- c(seSAMe1[[2]], seSAMe2[[2]],seSAMe3[[2]],seSAMe4[[2]],seSAMe5[[2]],seSAMe6[[2]],seSAMe7[[2]],seSAMe8[[2]],seSAMe9[[2]],seSAMe10[[2]])
con <- c(case1[[2]], case2[[2]], case3[[2]], case4[[2]], case5[[2]], case6[[2]], case7[[2]], case8[[2]], case9[[2]], case10[[2]])
mms <- c(mmstr1[[2]], mmstr2[[2]], mmstr3[[2]], mmstr4[[2]], mmstr5[[2]], mmstr6[[2]], mmstr7[[2]], mmstr8[[2]], mmstr9[[2]], mmstr10[[2]])
ses_df <- data.frame(cases = barcodes, SeSAMe = ses, Conumee = con, MethylMasteR = mms) %>% mutate(avg =  apply(ses_df[, c("SeSAMe", "Conumee", "MethylMasteR")], 1, mean)) %>% arrange(avg) 
plt_df <- ses_df %>% pivot_longer(names_to = "Technology", values_to = "Accuracy", cols = c("SeSAMe", "Conumee", "MethylMasteR")) %>% mutate(Technology = as.factor(Technology))

# Main plot.
p <- ggplot(plt_df, aes(x=Technology, y=Accuracy)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.05), size = 3, colour = '#01889F')

# Most accurate case and least accurate case plots.
leastAccCon <- matchingAlgo("TCGA-DX-A6BA", "Conumee", "LMS")
leastAccMethyl <- lmsProcessing("TCGA-DX-A6BA") %>% filter(!(type == "SNP"))
leastAccSesame <- softOutputProcessingLMS("TCGA-DX-A6BA") %>% filter(!(type == "SNP"))
leastAcc <- rbind(leastAccSesame, leastAccMethyl, leastAccCon) %>% dplyr::mutate(chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))
mostAccCon <- matchingAlgo("TCGA-DX-A7EN", "Conumee", "LMS")
mostAccMethyl <- lmsProcessing("TCGA-DX-A7EN") %>% filter(!(type == "SNP"))
mostAccSesame <- softOutputProcessingLMS("TCGA-DX-A7EN") %>% filter(!(type == "SNP"))
mostAcc <- rbind(mostAccCon, mostAccMethyl, mostAccSesame) %>% dplyr::mutate(chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))


mAcc <- plot_cnv_segments(mostAcc)

lAcc <- plot_cnv_segments(leastAcc)

# Conumee - LM
barcodesCM <- paste0("./Outputs/Conumee/LMData/", list.files("./Outputs/Conumee/LMData/", pattern = "\\.csv$") ) 

case1 <- fpCheck(matchingAlgoLM(barcodesCM[1]))
case2 <- fpCheck(matchingAlgoLM(barcodesCM[2]))
case3 <- fpCheck(matchingAlgoLM(barcodesCM[3]))
case4 <- fpCheck(matchingAlgoLM(barcodesCM[4]))
case5 <- fpCheck(matchingAlgoLM(barcodesCM[5]))
case6 <- fpCheck(matchingAlgoLM(barcodesCM[6]))
case7 <- fpCheck(matchingAlgoLM(barcodesCM[7])) # Issue - Resolved (NA Checking)
case8 <- fpCheck(matchingAlgoLM(barcodesCM[8])) # Issue - Resolved (NA Checking)
case9 <- fpCheck(matchingAlgoLM(barcodesCM[9])) # Issue - Resolved (NA Checking)

# SeSAMe - LM
barcodesSMLM <- list.files("./Outputs/SeSAMe/LM/") 
barcodesSMLM <- paste0("./Outputs/SeSAMe/LM/", barcodesSMLM[!(grepl(pattern = "bins_", barcodesSMLM))])

seSAMe1 <- fpCheck(softOutputProcessingLM(barcodesSMLM[1]))
seSAMe2 <- fpCheck(softOutputProcessingLM(barcodesSMLM[2]))
seSAMe3 <- fpCheck(softOutputProcessingLM(barcodesSMLM[3]))
seSAMe4 <- fpCheck(softOutputProcessingLM(barcodesSMLM[4]))
seSAMe5 <- fpCheck(softOutputProcessingLM(barcodesSMLM[5]))
seSAMe6 <- fpCheck(softOutputProcessingLM(barcodesSMLM[6]))
seSAMe7 <- fpCheck(softOutputProcessingLM(barcodesSMLM[7]))
seSAMe8 <- fpCheck(softOutputProcessingLM(barcodesSMLM[8])) 
seSAMe9 <- fpCheck(softOutputProcessingLM(barcodesSMLM[9])) 

# MethylMasteR - LM
barcodesMMLM <- paste0("./Outputs/MethylMaster/LM/", list.files("./Outputs/MethylMaster/LM", pattern = "autocorrected_regions.csv", recursive = T))
mmstr1 <- fpCheck(lmProcessing(barcodesMMLM[1]))
mmstr2 <- fpCheck(lmProcessing(barcodesMMLM[2]))
mmstr3 <- fpCheck(lmProcessing(barcodesMMLM[3]))
mmstr4 <- fpCheck(lmProcessing(barcodesMMLM[4]))
mmstr5 <- fpCheck(lmProcessing(barcodesMMLM[5]))
mmstr6 <- fpCheck(lmProcessing(barcodesMMLM[6]))
mmstr7 <- fpCheck(lmProcessing(barcodesMMLM[7])) # Issue
mmstr8 <- fpCheck(lmProcessing(barcodesMMLM[8])) # Issue
mmstr9 <- fpCheck(lmProcessing(barcodesMMLM[9])) # Issue

# LM Consolidation.
casesList <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6", "Sample 7", "Sample 8", "Sample 9")
ses <- c(seSAMe1[[2]], seSAMe2[[2]],seSAMe3[[2]],seSAMe4[[2]],seSAMe5[[2]],seSAMe6[[2]],seSAMe7[[2]],seSAMe8[[2]],seSAMe9[[2]])
con <- c(case1[[2]], case2[[2]], case3[[2]], case4[[2]], case5[[2]], case6[[2]], case7[[2]], case8[[2]], case9[[2]])
mms <- c(mmstr1[[2]], mmstr2[[2]], mmstr3[[2]], mmstr4[[2]], mmstr5[[2]], mmstr6[[2]], mmstr7[[2]], mmstr8[[2]], mmstr9[[2]])
ses_df <- data.frame(cases = casesList, SeSAMe = ses, Conumee = con, MethylMasteR = mms) %>% mutate(avg =  apply(ses_df[, c("SeSAMe", "Conumee", "MethylMasteR")], 1, mean)) %>% arrange(avg)
plt_df <- ses_df %>% pivot_longer(names_to = "Technology", values_to = "Accuracy", cols = c("SeSAMe", "Conumee", "MethylMasteR")) %>% mutate(Technology = as.factor(Technology))

# Main plot.
p <- ggplot(plt_df, aes(x=Technology, y=Accuracy)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.05), size = 3, colour = '#01889F')

leastAccCon <- matchingAlgoLM(barcodesCM[2])
leastAccMethyl <- lmProcessing(barcodesMMLM[2]) %>% filter(!(type == "SNP"))
leastAccSesame <- softOutputProcessingLM(barcodesSMLM[2]) %>% filter(!(type == "SNP"))
leastAcc <- rbind(leastAccSesame, leastAccMethyl, leastAccCon) %>% dplyr::mutate(chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))
mostAccCon <- matchingAlgoLM(barcodesCM[9])
mostAccMethyl <- lmProcessing(barcodesMMLM[9]) %>% filter(!(type == "SNP"))
mostAccSesame <- softOutputProcessingLM(barcodesSMLM[9]) %>% filter(!(type == "SNP"))
mostAcc <- rbind(mostAccCon, mostAccMethyl, mostAccSesame) %>% dplyr::mutate(chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))

mAcc <- plot_cnv_segments(mostAcc)

lAcc <- plot_cnv_segments(leastAcc)

