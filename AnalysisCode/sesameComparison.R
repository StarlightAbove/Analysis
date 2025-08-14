library(tidyverse)
softOutputProcessingLM <- function(outputDir){
  IDATSampleSheet <- read.csv("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
  chasFile <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_ChAS_CNVs.xlsx") %>% select(File, `Mean Log2Ratio`, Chromosome, Type, `Full Location`)
  fasst2file <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CNVs_LM/LM_FASST2_CNVs.xlsx")
  
  binSheet <- read.csv(str_replace(outputDir, "segments", replacement = "bins"))
  colnames(binSheet) <- c("loc", "value")
  binSheet <- tidyr::separate_wider_delim(binSheet, loc, names = c("chrom", "CpG"), delim = "-") %>% mutate(CpG = as.numeric(CpG)) %>% dplyr::mutate(
    chrom = str_remove_all(chrom, "chr")) %>% filter(!(chrom == "X") & !(chrom == "Y"))
  sesameOutput <- read.csv(outputDir) %>% dplyr::select(c("chrom", "loc.start", 
                  "loc.end", "seg.mean")) %>% dplyr::mutate(
                    chrom = str_remove_all(chrom, "chr")) %>% filter(
                      !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                        chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                          CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                          seg.mean < -0.3 ~ "Deletion",
                          TRUE ~ "Normal"), type = "SeSAMe") 
  # Need to filter it into chrom, loc.start, loc.end, seg.mean, CNVStatus & type
  caseNameCorrelation <- str_split(outputDir, "_")[[1]]
  sentrixID <- caseNameCorrelation[2]
  sentrixPos <- gsub(".csv", replacement = "", x = caseNameCorrelation[3])
  SUHMatch <- IDATSampleSheet %>% filter(Sentrix_ID == sentrixID, 
                                         Sentrix_Position == sentrixPos) %>% select(SUH) %>% dplyr::pull(SUH)
  
  chasFiltered <- chasFile[str_detect(chasFile$File, SUHMatch), ] %>% 
    select(-c(File)) %>% rename(chrom = Chromosome, CNVStatus = Type, seg.mean = `Mean Log2Ratio`) %>% 
    mutate(
      type = "SNP", 
      CNVStatus = case_when(CNVStatus == "Loss" ~ "Deletion", 
                          CNVStatus == "Gain" ~ "Amplification", 
                          TRUE ~ "Normal")
  )
  chasFiltered <-  tidyr::separate_wider_delim(tidyr::separate_wider_delim(chasFiltered, `Full Location`, names = c("chromosome", "loc"), delim = ":"), `loc`,
                    names = c("loc.start", "loc.end"), delim = "-") %>% select(-c("chromosome"))
  fasst2filtered <- fasst2file[str_detect(fasst2file$Sample, SUHMatch), ] %>%
                    select(c(`Probe Median`, `Event`, `Chromosome Region`)) %>% 
                    rename(seg.mean = `Probe Median`, CNVStatus = `Event`) %>%
                    mutate(CNVStatus = case_when(
                      CNVStatus == "CN Loss" ~ "Deletion", 
                      CNVStatus == "CN Gain" ~ "Amplification", 
                      TRUE ~ "Normal"
                    ))
  fasst2filtered <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(fasst2filtered, `Chromosome Region`, names = c("Chromosome", "loc"), delim = ":"), `loc`,
                                                names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") 
  fasst2filtered$loc.start <- as.numeric(str_replace_all(fasst2filtered$loc.start, ",", ""))
  fasst2filtered$loc.end <- as.numeric(str_replace_all(fasst2filtered$loc.end, ",", ""))
  fasst2filtered$chrom <- as.numeric(str_replace_all(fasst2filtered$chrom, "chr", ""))
  processing <- rbind(sesameOutput, chasFiltered, fasst2filtered) %>% 
    filter(!(CNVStatus == "Normal")) # Because of how ChAS and FASST2 work, a different paradigm must be used.
  processing$chrom <- as.numeric(processing$chrom)
  processing$loc.start <- as.numeric(processing$loc.start)
  processing$loc.end <- as.numeric(processing$loc.end)
  processing$seg.mean <- as.numeric(processing$seg.mean)
  
  # Binning into chromosomes.
  bins <- list()
  chromsomes <- unique(processing$chrom)
  for(crs in chromsomes){
    if(nrow(processing %>% dplyr::filter(chrom == crs)) > 0){
      bins[crs] <- list(processing %>% dplyr::filter(chrom == crs))
    } 
  }
  Filter(Negate(is.null), bins)
}


# Statistical function.
stats <- function(binnedList){
  mets <- data.frame(chroms = 1:22, acc = c(0), fp = c(0), fn = c(0), exists = F)
  
  for(bin in binnedList){
    mets[bin$chrom, ]$exists <- T
  }
  
  mets <- mets %>% filter(exists == T) %>% select(-c("exists"))
  
  for(bin in binnedList){
    segs <- c()
    softs <- bin %>% filter(!(type == "SNP"))
    
    if(nrow(softs) == 0){
      mets[unique(mets$chrom) == unique(bin$chrom), ]$fn <- mets[unique(mets$chrom) == unique(bin$chrom), ]$fn + nrow(bin %>% filter(type == "SNP"))
    }
    
    if(nrow(softs) > 0 & nrow(bin %>% filter(type == "SNP")) == 0){
      mets[unique(mets$chrom) == unique(bin$chrom), ]$fp <- mets[unique(mets$chrom) == unique(bin$chrom), ]$fp + nrow(bin %>% filter(!(type == "SNP")))
    }
    
    if(nrow(softs) > 0 & nrow(bin %>% filter(type == "SNP")) > 0){
      for(sf in 1:nrow(softs)){
        
        seg <- softs[sf, ]
        segArr <- seg$loc.end - seg$loc.start
        
        # Sub-setting SNPs.
        snpSubset <- bin %>% filter(type == "SNP") %>% filter(loc.start > seg$loc.start & loc.end < seg$loc.end)
        topEdge <- bin %>% filter(loc.start < seg$loc.start & loc.end > seg$loc.start & ((loc.end - seg$loc.start)/(loc.end - loc.start) > 0.5), type == "SNP")
        if(nrow(topEdge) >= 1){
          topEdge$loc.start <- seg$loc.start # To normalize it to within the segment.
          if(nrow(topEdge) >= 1){
            topEdge <- topEdge %>% arrange(desc(loc.end))
            topEdge <- topEdge[1, ]
          }
        }
        
        bottomEdge <- bin %>% filter(loc.start < seg$loc.end & loc.end > seg$loc.end & loc.start > seg$loc.start & ((seg$loc.end - loc.start)/(loc.end - loc.start) > 0.5), type == "SNP")
        if(nrow(bottomEdge) >= 1){
          bottomEdge$loc.end <- seg$loc.end # Normalize it to within the segment.
          if(nrow(bottomEdge) >= 1){
            bottomEdge <- bottomEdge %>% arrange(desc(loc.end))
            bottomEdge <- bottomEdge[1, ]
          }
        }
        
        if(nrow(bottomEdge) >= 1 | nrow(topEdge) >= 1){
          snpSubset <- rbind(snpSubset, topEdge, bottomEdge) %>% dplyr::arrange(loc.start)
        }
        
        
        if(nrow(snpSubset) == 1){
          geneLen <- snpSubset[1, ]$loc.end - snpSubset[1, ]$loc.start
          segs <- c(segs, geneLen/segArr)
        }
        
        if(nrow(snpSubset) == 0){
          segs <- c(segs, 0)
        }
        
        if(nrow(snpSubset) > 1){
          deletion <- c()
          amplification <- c()
          i <- 1
          while(i < nrow(snpSubset)){
            startBlock <- snpSubset[i, ]
            j <- i + 1
            endBlock <- snpSubset[j, ]
            while(j <= nrow(snpSubset)){
              endBlock <- snpSubset[j, ]
              if(!(startBlock$CNVStatus == endBlock$CNVStatus)){
                j <- j - 1
                if(startBlock$CNVStatus == "Deletion"){
                  deletion <- c(deletion, endBlock$loc.end - startBlock$loc.start)
                }
                if(startBlock$CNVStatus == "Amplification"){
                  amplification <- c(amplification, endBlock$loc.end - startBlock$loc.start)
                }
                break
              }
              j <- j + 1
            }
          }
          
          if(seg$CNVStatus == "Deletion"){
            segs <- c(segs, sum(deletion)/segArr)
          }
          if(seg$CNVStatus == "Amplification"){
            segs <- c(segs, sum(amplification)/segArr)
          }
        }
        
      }
      mets[unique(mets$chrom) == unique(bin$chrom), ]$acc <- sum(segs)/length(segs)
    }
    
  }
  
  list(byChromosome = mets, caseAcc = (sum(t$acc) + sum(22 - nrow(t)))/22)
  
}
  
barcodes <- list.files("./Outputs/SeSAMe/LM") 
barcodes <- paste0("./Outputs/SeSAMe/LM/", barcodes[!(grepl(pattern = "bins_", barcodes))])

t1 <- softOutputProcessingLM(barcodes[1])
t1stat <- stats(t1) # Clear.
t2 <- softOutputProcessingLM(barcodes[2])
t2stat <- stats(t2)
t3 <- softOutputProcessingLM(barcodes[3])
t3stat <- stats(t3) # Not clear.
t4 <- softOutputProcessingLM(barcodes[4])
t4stat <- stats(t4)
t5 <- softOutputProcessingLM(barcodes[5])
t5stat <- stats(t5)
t6 <- softOutputProcessingLM(barcodes[6])
t6stat <- stats(t6)
t7 <- softOutputProcessingLM(barcodes[7])
t7stat <- stats(t7) # Not clear.
t8 <- softOutputProcessingLM(barcodes[8])
t8stat <- stats(t8)
t9 <- softOutputProcessingLM(barcodes[9])
t9stat <- stats(t9)



