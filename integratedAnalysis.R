# Document to analyse datasets on different software, as previous file was more experimental.

# HM450 controls.
library(GEOquery)
library(tidyverse)
library(R.utils)

Ctrls <- c("GSM3015205", "GSM3015206", "GSM3015207", "GSM3015208","GSM3015209","GSM3015210")
setwd("~/Work/MethylVignettes/GDC Analysis/Controls/HM450Controls")
for(i in Ctrls) getGEOSuppFiles(GEO = i, makeDirectory = T)
rm(i)
setwd("~/Work/MethylVignettes/GDC Analysis")
for(ct in Ctrls){
  f <- list.files(path = paste0("./Controls/HM450Controls/", ct, collapse="/"), full.names = F)
  fldrName <- str_remove(f[1], "_Grn.idat.gz")
  file.rename(paste0("./Controls/HM450Controls/", ct, collapse="/"), fldrName)
  
}

f <- list.files(path = list.dirs(path = "./Controls/HM450Controls"), pattern = "\\.idat.gz$", full.names = TRUE)
for(fs in f){
  gunzip(fs,remove = F)
}
controls <- list.files(path = list.dirs(path = "./Controls/HM450Controls"), pattern = "\\.idat$", full.names = TRUE)
dest <- "./Controls/HM450Controls/ControlFile"
for(ctls in controls){
  file.copy(ctls, dest)
}

# Leiomyoma data -> LM_SNP_EPIC_array_data
# Conumee - Processed.
library(conumee2)
RGset <- read.metharray.exp("./LabData/LM_SNP_EPIC_array_data", recursive = T, verbose = T)
MSet <- preprocessIllumina(RGset)
CSet <- read.metharray.exp("./Controls", recursive = T, verbose = T)
CSet <- preprocessIllumina(CSet)
anno <- CNV.create_anno(array_type = c("450k", "EPIC"))
load.data <- CNV.load(MSet)
load.controls <- CNV.load(CSet)
x <- CNV.fit(load.data, ref = load.controls, anno)
x <- CNV.bin(x)
# x <- CNV.detail(x) (Can be added if there is certain focus on certain loci)
x <- CNV.segment(x)
segments <- CNV.write(x, what = "segments")
bins <- CNV.write(x, what = "bins")
write.csv(bins, "./Outputs/Conumee/LMData/bins.csv")
for(segment in segments){
  write.csv(segment, paste0("./Outputs/Conumee/LMData/", segment[["ID"]][1], ".csv"))
}

# SeSAMe - Processed.
library(sesame)
library(sesameData)
library(reshape2)
idat_dir <- "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM"
sdfs <- openSesame(idat_dir, func = NULL)

i <- 1
for(sdf in sdfs){
  segs <- sdf
  segments <- cnSegmentation(segs)
  segmentalSignals <- segments[["seg.signals"]]
  bins <- segments[["bin.signals"]]
  bins <- reshape2::melt(bins)
  write.csv(bins, paste0("./Outputs/SeSAMe/LM/bins_",names(sdfs[i]),".csv"), row.names = TRUE)
  write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/LM/segments_", 
                                     names(sdfs[i]), ".csv"), row.names = FALSE)
  i <- i + 1
}

# Epicopy.
library(Epicopy)
# Epicopy seems to not be able to process EPIC & EPICv2 files. Therefore, this data is skipped. 

# ChAMP.
library(ChAMP)

# Converting xlsx file to csv.
xlsx <- readxl::read_excel("./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/design_9LM_EPIC_array.xlsx")
csv <- write.csv(xlsx, "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv")
# ChAMP cannot handle EPICv2.
myLoad <- champ.load(directory = "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files", arraytype = "EPIC", method = "ChAMP")
myCNA <- champ.CNA(control = F,arraytype = "EPICv1", groupFreqPlots=FALSE)

# MethylMastR. - Processed.
library(MethylMasteR)
library(sesameData)
sesameDataCacheAll()
library(tidyverse)
input.dir <- paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files")
output.dir <- paste0(getwd(),"/Outputs/MethylMaster/LM")
routine.run <- "custom"
EPICSheet <- read.csv(paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv"))
sampSheet <- data.frame(Sample_Name = paste0(EPICSheet$Sentrix_ID, "_", EPICSheet$Sentrix_Position),
                        Sample_Group = "tumor",
                        Sample_Plate = EPICSheet$Sample_Plate,
                        Pool_ID = NA,
                        Sentrix_Position = EPICSheet$Sentrix_Position,
                        Sentrix_ID = EPICSheet$Sentrix_ID,
                        Basename = paste0(EPICSheet$Sentrix_ID, "_", EPICSheet$Sentrix_Position),
                        Platform = "EPIC",
                        Batch = 1,
                        Tumor = 75.76724742)
write_csv(sampSheet, "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/MethylSheet.csv") 
sample.sheet.path <- paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/MethylSheet.csv")

methyl_master(
  routine                   = routine.run, #The routine to run
  input.dir                 = input.dir, #The input (idat.files) directory
  output.dir                = output.dir, #The output directory
  sample.sheet.path         = sample.sheet.path, #The path to the MethylMasteR sample sheet
  r.lib.path                = .libPaths()[1], #The path to the R Library path
  file.sep                  = "\\\\", #For windows or "/" for Linux
  create.dir                = TRUE, #Whether to cretae directory if does not 
  #exist?
  save.seg                  = TRUE, #Whether to save segmentation results
  n.cores                   = 1, #Multicore does not work for all routines 
  #on all operating systems
  os.type                   = "linux", #Or "linux"
  proj                      = "TCGA-KIRC", #"TCGA-BLCA" etc.
  visualize                 = TRUE, #Whether to output plots,
  visualize.individual      = FALSE, #Whether to output indvidual sample plots,
  #only works for routine sesame
  reference                 = "internal", #"comparison" or 'internal"
  reference.name            = NA, #For Epicopy use NA for median, 
  #"all" is not currently supported 
  #by Epicopy
  comparison                = c("tumor","normal"), #Always required treatment 
  #first then contro can also be more specific when 
  #designing sample sheet and use values like
  # c("tumor_male","cord_male etc") etc.
  #Note: in routines where internal reference is 
  #used, second argument is ignored
  form.thresholds           = NULL, #Used to calculate final CNV state. 
  #If NULL, equation is is used;
  #otherwise, specify threshold vector 
  #of lower and upper beta values such
  #as c(-0.3,0.3), we used c(-0.2,0.2) 
  #in the paper
  overlap.density           = 0.1, #For combining final CNV calls for confidence
  sesame.data.cache         = "EPIC", #The default sesame reference platform
  sesame.data.normal        = 'EPIC.5.normal', #The default sesame and hm450 
  #internal reference samples
  genome.version            = "hg19", #Or can set to "hg38"
  hm450.workflow            = "B", #The HM450 subworkflow to use - only B 
  #is running currently.
  #"A" no correction,
  #"B" median correction (default),
  #"C" run Conumee
  champ.padj                = 0.05, #padj to filter champ results
  champ.control             = FALSE, #run champ.control etc.
  champ.run.combat          = FALSE, #run champ.run.combat etc.
  champ.run.dmp             = FALSE, #If only one pheno var must = FALSE
  champ.run.dmr             = FALSE, #If only one pheno var must = FALSE
  champ.run.block           = FALSE, #If only one pheno var must = FALSE
  champ.run.gsea            = FALSE, #Requires dmp and dmr results
  champ.run.epimod          = FALSE, #If only one pheno var must = FALSE
  epi.run.gistic            = TRUE, #Whether to Run GISTIC in Epicopy workflow
  olaps.split.field         = "Sample_ID", #Split field to ise during overlaps
  #Don't change unless you know what 
  #you are doing
  estimate.recurrence       = TRUE, #Estimate recursion to produce p values when 
  #finding overlaps with population_ranges 
  #functions
  ov.pvalue                 = 0.05, #pvalue threshold for overlaps identified
  ov.keep.extra.columns     = TRUE, #Keep extra metadata columns when finding 
  #overlaps
  simplify.reduce           = weightedmean #Equation to use during reduction
)




# Leiomyosarcoma.
# Conumee. - Processed.
f <- unique(str_remove(list.files("./Controls/", recursive = T, pattern = "*.idat$"), pattern = "_Grn.idat|_Red.idat"))
RGset <- read.metharray.exp("./LabData", recursive = T, verbose = T)
MSet <- preprocessIllumina(RGset)
CSet <- read.metharray.exp("./Controls", recursive = T, verbose = T)
CSet <- preprocessIllumina(CSet)
anno <- CNV.create_anno(array_type = c("450k", "EPIC"))
load.data <- CNV.load(MSet)
load.controls <- CNV.load(CSet)
x <- CNV.fit(load.data, ref = load.controls, anno)
x <- CNV.bin(x)
# x <- CNV.detail(x) (If detailed analysis wanted to be done.)
x <- CNV.segment(x)
segments <- CNV.write(x, what = "segments")
bins <- CNV.write(x, what = "bins")
write.csv(bins, "./Outputs/Conumee/PrzybylLabData/bins.csv")
for(segment in segments){
  write.csv(segment, paste0("./Outputs/Conumee/PrzybylLabData/", segment[["ID"]][1], ".csv"))
}

# SeSAMe. - Processed.
library(sesame)
library(sesameData)
library(reshape2)
idat_dir <- "./cases/SingleFileCase"
sdfs <- openSesame(idat_dir, func = NULL)
Control <- openSesame("./Controls/HM450Controls/ControlFile", func = NULL)


i <- 1
for(sdf in sdfs){
  segs <- sdf
  segments <- cnSegmentation(segs, sdfs.normal = Control)
  segmentalSignals <- segments[["seg.signals"]]
  bins <- segments[["bin.signals"]]
  bins <- reshape2::melt(bins)
  write.csv(bins, paste0("./Outputs/SeSAMe/LMS/bins_",names(sdfs[i]),".csv"), row.names = TRUE)
  write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/LMS/segments_", 
                                     names(sdfs[i]), ".csv"), row.names = FALSE)
  i <- i + 1
}

# ChAMP. - Processed but needs to be re-run to get actual textual data out of it.
library(ChAMP)
# Building a sample sheet & files according to Illumina requirements.
TCGASheet <- read_tsv(
  "./cases/SingleFileCase/Sample_Sheet.tsv"
) %>% dplyr::select(c("File Name", "Case ID"))
SampleSheet <- read_csv("./cases/SingleFileCase/SampleSheet.csv")
for(i in 1:nrow(SampleSheet)){
  selection <- SampleSheet[i, ]
  renamer <- TCGASheet %>% filter(`Case ID` == selection$Sample_Name)
  for(j in 1:nrow(renamer)){
    newName <- paste0(selection$Sentrix_ID, "_", selection$Sentrix_Position, "_", str_split(renamer[j, ]$`File Name`, pattern = "_")[[1]][3])
    print(newName)
    file.rename(from = paste0("./cases/SingleFileCase/", renamer[j, ]$`File Name`), to = paste0("./cases/SingleFileCase/", newName))
  }
}



myLoad <- champ.load(directory = "./cases/SingleFileCase", arraytype = "HM450", method = "ChAMP")
myCNA <- champ.CNA(control = F,arraytype = "HM450", groupFreqPlots=T)

# Epicopy. This package only works with HM450K. 
library(Epicopy)
library(minfi)
epi_ss <- read.metharray.sheet("./cases/SingleFileCase", pattern = "SampleSheet.csv")
epi_rg <- read.metharray.exp(targets = epi_ss)
epi_lrr <- getLRR(rgSet = epi_rg, Normals = NA)

# MethylMasteR. - Processed.
library(MethylMasteR)
library(sesameData)
sesameDataCacheAll()
library(tidyverse)
input.dir <- paste0(getwd(), "/cases/SingleFileCase")
output.dir <- paste0(getwd(),"/Outputs/MethylMaster/LMS")
# sampSheet <- data.frame(Sample_Name = paste0(SampleSheet$Sentrix_ID, "_", SampleSheet$Sentrix_Position), 
# Sample_Group = SampleSheet$Sample_Group,
# Sample_Plate = NA,
# Pool_ID = NA,
# Sentrix_Position = SampleSheet$Sentrix_Position,
# Sentrix_ID = SampleSheet$Sentrix_ID,
# Basename = paste0(SampleSheet$Sentrix_ID, "_", SampleSheet$Sentrix_Position),
# Platform = "HM450",
# Batch = 1,
# Tumor = 75.76724742) # Building a MethylMaster compliant datasheet.
sampSheet <- sampSheet %>% mutate(caseNames = SampleSheet$Sample_Name)
write_csv(sampSheet, "./cases/SingleFileCase/SampleSheet.csv")
sample.sheet.path <- paste0(getwd(), "/cases/SingleFileCase/SampleSheet.csv")
routine.run <- "custom"

methyl_master(
  routine                   = routine.run, #The routine to run
  input.dir                 = input.dir, #The input (idat.files) directory
  output.dir                = output.dir, #The output directory
  sample.sheet.path         = sample.sheet.path, #The path to the MethylMasteR sample sheet
  r.lib.path                = .libPaths()[1], #The path to the R Library path
  file.sep                  = "\\\\", #For windows or "/" for Linux
  create.dir                = TRUE, #Whether to cretae directory if does not 
  #exist?
  save.seg                  = TRUE, #Whether to save segmentation results
  n.cores                   = 1, #Multicore does not work for all routines 
  #on all operating systems
  os.type                   = "linux", #Or "linux"
  proj                      = "TCGA-KIRC", #"TCGA-BLCA" etc.
  visualize                 = TRUE, #Whether to output plots,
  visualize.individual      = FALSE, #Whether to output indvidual sample plots,
  #only works for routine sesame
  reference                 = "internal", #"comparison" or 'internal"
  reference.name            = NA, #For Epicopy use NA for median, 
  #"all" is not currently supported 
  #by Epicopy
  comparison                = c("tumor","normal"), #Always required treatment 
  #first then contro can also be more specific when 
  #designing sample sheet and use values like
  # c("tumor_male","cord_male etc") etc.
  #Note: in routines where internal reference is 
  #used, second argument is ignored
  form.thresholds           = NULL, #Used to calculate final CNV state. 
  #If NULL, equation is is used;
  #otherwise, specify threshold vector 
  #of lower and upper beta values such
  #as c(-0.3,0.3), we used c(-0.2,0.2) 
  #in the paper
  overlap.density           = 0.1, #For combining final CNV calls for confidence
  sesame.data.cache         = "EPIC", #The default sesame reference platform
  sesame.data.normal        = 'EPIC.5.normal', #The default sesame and hm450 
  #internal reference samples
  genome.version            = "hg19", #Or can set to "hg38"
  hm450.workflow            = "B", #The HM450 subworkflow to use - only B 
  #is running currently.
  #"A" no correction,
  #"B" median correction (default),
  #"C" run Conumee
  champ.padj                = 0.05, #padj to filter champ results
  champ.control             = FALSE, #run champ.control etc.
  champ.run.combat          = FALSE, #run champ.run.combat etc.
  champ.run.dmp             = FALSE, #If only one pheno var must = FALSE
  champ.run.dmr             = FALSE, #If only one pheno var must = FALSE
  champ.run.block           = FALSE, #If only one pheno var must = FALSE
  champ.run.gsea            = FALSE, #Requires dmp and dmr results
  champ.run.epimod          = FALSE, #If only one pheno var must = FALSE
  epi.run.gistic            = TRUE, #Whether to Run GISTIC in Epicopy workflow
  olaps.split.field         = "Sample_ID", #Split field to ise during overlaps
  #Don't change unless you know what 
  #you are doing
  estimate.recurrence       = TRUE, #Estimate recursion to produce p values when 
  #finding overlaps with population_ranges 
  #functions
  ov.pvalue                 = 0.05, #pvalue threshold for overlaps identified
  ov.keep.extra.columns     = TRUE, #Keep extra metadata columns when finding 
  #overlaps
  simplify.reduce           = weightedmean #Equation to use during reduction
)

