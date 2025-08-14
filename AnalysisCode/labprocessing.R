library(tidyverse)
library(R.utils)

# Leiomyosarcoma data -> LMS_SNP_EPIC_array_data
# Conumee. - Processed.
library(conumee2)

RGset <- read.metharray.exp("./LabData/LMS_SNP_EPIC_array_data", recursive = T, verbose = T)
MSet <- preprocessIllumina(RGset)
CSet <- read.metharray.exp("./LabData/Normal_smooth_muscle_EPIC_data", recursive = T, verbose = T)
CSet <- preprocessIllumina(CSet)
anno <- CNV.create_anno(array_type = "EPIC")
load.data <- CNV.load(MSet)
load.controls <- CNV.load(CSet)
x <- CNV.fit(load.data, ref = load.controls, anno)
x <- CNV.bin(x)
# x <- CNV.detail(x) (Can be added if there is certain focus on certain loci)
x <- CNV.segment(x)
segments <- CNV.write(x, what = "segments")
bins <- CNV.write(x, what = "bins")

write.csv(bins, "./Outputs/Conumee/LabLMS/bins.csv")
for(segment in segments){
  write.csv(segment, paste0("./Outputs/Conumee/LabLMS/", segment[["ID"]][1], ".csv"))
}

# SeSAMe.
library(sesame)
library(sesameData)
library(reshape2)
idat_dir <- "./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files"
sdfs <- openSesame(idat_dir, func = NULL)

i <- 1
for(sdf in sdfs){
  segs <- sdf
  segments <- cnSegmentation(segs)
  segmentalSignals <- segments[["seg.signals"]]
  bins <- segments[["bin.signals"]]
  bins <- reshape2::melt(bins)
  write.csv(bins, paste0("./Outputs/SeSAMe/LabLMS/bins_",names(sdfs[i]),".csv"), row.names = TRUE)
  write.csv(segmentalSignals, paste0("./Outputs/SeSAMe/LabLMS/segments_", 
                                     names(sdfs[i]), ".csv"), row.names = FALSE)
  i <- i + 1
}


# MethylMasteR. - Finish in Rocker session.
library(MethylMasteR)
library(sesameData)
sesameDataCacheAll()
library(tidyverse)
input.dir <- paste0(getwd(), "/LabData/MethylMasteRCombined")
output.dir <- paste0(getwd(),"/Outputs/MethylMaster/LabLMS")
routine.run <- "custom"
sample.sheet.path <- paste0(getwd(), "/LabData/MethylMasteRCombined/Sample_Sheet_Combined.csv")

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
  form.thresholds           = c(-0.2, 0.2), #Used to calculate final CNV state. 
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

