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

# Lab LMS Consolidation.
# Conumee - L:LMS
correlative <- read_csv("LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
case1 <- fpCheck(labLMSProc(correlative[1], "Conumee"))
case2 <- fpCheck(labLMSProc(correlative[2], "Conumee"))
case3 <- fpCheck(labLMSProc(correlative[3], "Conumee"))
case4 <- fpCheck(labLMSProc(correlative[4], "Conumee"))
case5 <- fpCheck(labLMSProc(correlative[5], "Conumee"))
case6 <- fpCheck(labLMSProc(correlative[6], "Conumee"))
case7 <- fpCheck(labLMSProc(correlative[7], "Conumee")) # Issue - Resolved (NA Checking)
case8 <- fpCheck(labLMSProc(correlative[8], "Conumee")) # Issue - Resolved (NA Checking)
case9 <- fpCheck(labLMSProc(correlative[9], "Conumee"))
case10 <- fpCheck(labLMSProc(correlative[10], "Conumee"))
case11 <- fpCheck(labLMSProc(correlative[11], "Conumee"))
case12 <- fpCheck(labLMSProc(correlative[12], "Conumee"))
case13 <- fpCheck(labLMSProc(correlative[13], "Conumee"))

# Sesame - L:LMS
correlative <- read_csv("LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
sm1 <- fpCheck(labLMSProc(correlative[1], "Sesame"))
sm2 <- fpCheck(labLMSProc(correlative[2], "Sesame"))
sm3 <- fpCheck(labLMSProc(correlative[3], "Sesame"))
sm4 <- fpCheck(labLMSProc(correlative[4], "Sesame"))
sm5 <- fpCheck(labLMSProc(correlative[5], "Sesame"))
sm6 <- fpCheck(labLMSProc(correlative[6], "Sesame"))
sm7 <- fpCheck(labLMSProc(correlative[7], "Sesame")) # Issue - Resolved (NA Checking)
sm8 <- fpCheck(labLMSProc(correlative[8], "Sesame")) # Issue - Resolved (NA Checking)
sm9 <- fpCheck(labLMSProc(correlative[9], "Sesame"))
sm10 <- fpCheck(labLMSProc(correlative[10], "Sesame"))
sm11 <- fpCheck(labLMSProc(correlative[11], "Sesame"))
sm12 <- fpCheck(labLMSProc(correlative[12], "Sesame"))
sm13 <- fpCheck(labLMSProc(correlative[13], "Sesame"))

# MethylMaster - L:LMS
correlative <- read_csv("LabData/LMS_SNP_EPIC_array_data/correlative.csv")$STT
mm1 <- fpCheck(labLMSProc(correlative[1], "MethylMaster"))
mm2 <- fpCheck(labLMSProc(correlative[2], "MethylMaster"))
mm3 <- fpCheck(labLMSProc(correlative[3], "MethylMaster"))
mm4 <- fpCheck(labLMSProc(correlative[4], "MethylMaster"))
mm5 <- fpCheck(labLMSProc(correlative[5], "MethylMaster"))
mm6 <- fpCheck(labLMSProc(correlative[6], "MethylMaster"))
mm7 <- fpCheck(labLMSProc(correlative[7], "MethylMaster")) # Issue - Resolved (NA Checking)
mm8 <- fpCheck(labLMSProc(correlative[8], "MethylMaster")) # Issue - Resolved (NA Checking)
mm9 <- fpCheck(labLMSProc(correlative[9], "MethylMaster"))
mm10 <- fpCheck(labLMSProc(correlative[10], "MethylMaster"))
mm11 <- fpCheck(labLMSProc(correlative[11], "MethylMaster"))
mm12 <- fpCheck(labLMSProc(correlative[12], "MethylMaster"))
mm13 <- fpCheck(labLMSProc(correlative[13], "MethylMaster"))

# L:LMS Consolidation.
casesList <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", 
               "Sample 6", "Sample 7", "Sample 8", "Sample 9", "Sample 10", 
               "Sample 11", "Sample 12", "Sample 13")
ses <- c(sm1[[2]], sm2[[2]],sm3[[2]],sm4[[2]],sm5[[2]],sm6[[2]],sm7[[2]],sm8[[2]],sm9[[2]], sm10[[2]], sm11[[2]], sm12[[2]], sm13[[2]])
con <- c(case1[[2]], case2[[2]], case3[[2]], case4[[2]], case5[[2]], case6[[2]], case7[[2]], case8[[2]], case9[[2]], case10[[2]], case11[[2]], case12[[2]], case13[[2]])
mms <- c(mm1[[2]], mm2[[2]],mm3[[2]],mm4[[2]],mm5[[2]],mm6[[2]],mm7[[2]],mm8[[2]],mm9[[2]], mm10[[2]], mm11[[2]], mm12[[2]], mm13[[2]])
ses_df <- data.frame(cases = casesList, SeSAMe = ses, Conumee = con, MethylMasteR = mms)
ses_df <- ses_df %>% dplyr::mutate(avg =  apply(ses_df[, c("SeSAMe", "Conumee", "MethylMasteR")], 1, mean)) %>% arrange(avg)
plt_df <- ses_df %>% pivot_longer(names_to = "Technology", values_to = "Accuracy", cols = c("SeSAMe", "Conumee", "MethylMasteR")) %>% mutate(Technology = as.factor(Technology))
p <- ggplot(plt_df, aes(x=Technology, y=Accuracy)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.05), size = 3, colour = '#01889F')

leastAccCon <- labLMSProc(correlative[8], "Conumee") %>% filter(!(type == "SNP"))
leastAccMethyl <- labLMSProc(correlative[8], "MethylMaster") %>% filter(!(type == "SNP"))
leastAccSesame <- labLMSProc(correlative[8], "Sesame") %>% filter(!(type == "SNP"))
leastAcc <- rbind(leastAccSesame, leastAccMethyl, leastAccCon) %>% dplyr::mutate(chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))
mostAccCon <- labLMSProc(correlative[2], "Conumee") %>% filter(!(type == "SNP"))
mostAccMethyl <- labLMSProc(correlative[2], "MethylMaster") %>% filter(!(type == "SNP"))
mostAccSesame <- labLMSProc(correlative[2], "Sesame") %>% filter(!(type == "SNP"))
mostAcc <- rbind(mostAccCon, mostAccMethyl, mostAccSesame) %>% dplyr::mutate(chrom = as.numeric(chrom)) %>% filter(!(is.na(chrom)))

mAcc <- plot_cnv_segments(mostAcc)

lAcc <- plot_cnv_segments(leastAcc)



# All Manhattan Plots.
# LMS Plots.
barcodes <- as.list(read.csv("./cases/SingleFileCase/SampleSheet.csv") %>% select("caseNames"))[[1]]
smpLMS1 <- rbind(matchingAlgo(barcodes[1], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[1]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[1]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS1)
smpLMS2 <- rbind(matchingAlgo(barcodes[2], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[2]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[2]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS2)
smpLMS3 <- rbind(matchingAlgo(barcodes[3], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[3]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[3]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS3)
smpLMS4 <- rbind(matchingAlgo(barcodes[4], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[4]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[4]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS4)
smpLMS5 <- rbind(matchingAlgo(barcodes[5], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[5]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[5]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS5)
smpLMS6 <- rbind(matchingAlgo(barcodes[6], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[6]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[6]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS6)
smpLMS7 <- rbind(matchingAlgo(barcodes[7], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[7]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[7]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS7)
smpLMS8 <- rbind(matchingAlgo(barcodes[8], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[8]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[8]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS8)
smpLMS9 <- rbind(matchingAlgo(barcodes[9], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[9]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[9]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS9)
smpLMS10 <- rbind(matchingAlgo(barcodes[10], "Conumee", "LMS"), (softOutputProcessingLMS(barcodes[10]) %>% filter(!(type == "SNP"))), (lmsProcessing(barcodes[10]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLMS10)

# LM Plots.
barcodesCM <- paste0("./Outputs/Conumee/LMData/", list.files("./Outputs/Conumee/LMData/", pattern = "\\.csv$") ) 
barcodesSMLM <- list.files("./Outputs/SeSAMe/LM/") 
barcodesSMLM <- paste0("./Outputs/SeSAMe/LM/", barcodesSMLM[!(grepl(pattern = "bins_", barcodesSMLM))])
barcodesMMLM <- paste0("./Outputs/MethylMaster/LM/", list.files("./Outputs/MethylMaster/LM", pattern = "autocorrected_regions.csv", recursive = T))

smpLM1 <- rbind(matchingAlgoLM(barcodesCM[1]), (softOutputProcessingLM(barcodesSMLM[1]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[1]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM1)
smpLM2 <- rbind(matchingAlgoLM(barcodesCM[2]), (softOutputProcessingLM(barcodesSMLM[2]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[2]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM2)
smpLM3 <- rbind(matchingAlgoLM(barcodesCM[3]), (softOutputProcessingLM(barcodesSMLM[3]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[3]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM3)
smpLM4 <- rbind(matchingAlgoLM(barcodesCM[4]), (softOutputProcessingLM(barcodesSMLM[4]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[4]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM4)
smpLM5 <- rbind(matchingAlgoLM(barcodesCM[5]), (softOutputProcessingLM(barcodesSMLM[5]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[5]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM5)
smpLM6 <- rbind(matchingAlgoLM(barcodesCM[6]), (softOutputProcessingLM(barcodesSMLM[6]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[6]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM6)
smpLM7 <- rbind(matchingAlgoLM(barcodesCM[7]), (softOutputProcessingLM(barcodesSMLM[7]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[7]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM7)
smpLM8 <- rbind(matchingAlgoLM(barcodesCM[8]), (softOutputProcessingLM(barcodesSMLM[8]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[8]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM8)
smpLM9 <- rbind(matchingAlgoLM(barcodesCM[9]), (softOutputProcessingLM(barcodesSMLM[9]) %>% filter(!(type == "SNP"))), (lmProcessing(barcodesMMLM[9]) %>% filter(!(type == "SNP"))))
plot_cnv_segments(smpLM9)


# Lab LMS Plots.
smp1 <- rbind(labLMSProc(correlative[1], "Sesame"), rbind(labLMSProc(correlative[1], "Conumee"), labLMSProc(correlative[1], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp1)
smp2 <- rbind(labLMSProc(correlative[2], "Sesame"), rbind(labLMSProc(correlative[2], "Conumee"), labLMSProc(correlative[2], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp2)
smp3 <- rbind(labLMSProc(correlative[3], "Sesame"), rbind(labLMSProc(correlative[3], "Conumee"), labLMSProc(correlative[3], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp3)
smp4 <- rbind(labLMSProc(correlative[4], "Sesame"), rbind(labLMSProc(correlative[4], "Conumee"), labLMSProc(correlative[4], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp4)
smp5 <- rbind(labLMSProc(correlative[5], "Sesame"), rbind(labLMSProc(correlative[5], "Conumee"), labLMSProc(correlative[5], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp5)
smp6 <- rbind(labLMSProc(correlative[6], "Sesame"), rbind(labLMSProc(correlative[6], "Conumee"), labLMSProc(correlative[6], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp6)
smp7 <- rbind(labLMSProc(correlative[7], "Sesame"), rbind(labLMSProc(correlative[7], "Conumee"), labLMSProc(correlative[7], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp7)
smp8 <- rbind(labLMSProc(correlative[8], "Sesame"), rbind(labLMSProc(correlative[8], "Conumee"), labLMSProc(correlative[8], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp8)
smp9 <- rbind(labLMSProc(correlative[9], "Sesame"), rbind(labLMSProc(correlative[9], "Conumee"), labLMSProc(correlative[9], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp9)
smp10 <- rbind(labLMSProc(correlative[10], "Sesame"), rbind(labLMSProc(correlative[10], "Conumee"), labLMSProc(correlative[10], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp10)
smp11 <- rbind(labLMSProc(correlative[11], "Sesame"), rbind(labLMSProc(correlative[11], "Conumee"), labLMSProc(correlative[11], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp11)
smp12 <- rbind(labLMSProc(correlative[12], "Sesame"), rbind(labLMSProc(correlative[12], "Conumee"), labLMSProc(correlative[12], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp12)
smp13 <- rbind(labLMSProc(correlative[13], "Sesame"), rbind(labLMSProc(correlative[13], "Conumee"), labLMSProc(correlative[13], "MethylMaster") %>% dplyr::filter(!(type == "SNP"))))
plot_cnv_segments(smp13)



