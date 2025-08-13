# Processing all methylation datapoints into CNVs.
library("conumee2")
library("sesame")
library("minfi")
library("readr")
library(MethylMasteR)
library("ChAMP")


# Control files from SERIES GSE192353.
dir.create("Controls")
Ctrls <- c("GSM5745337", "GSM5745338","GSM5745339",	"GSM5745340",	"GSM5745341",	"GSM5745342",	"GSM5745343",	
           "GSM5745344",	"GSM5745345",	"GSM5745346",	"GSM5745347",	"GSM5745348",	
           "GSM5745349",	"GSM5745350",	"GSM5745351",	"GSM5745352")
# list.files(path = "./Controls", pattern = "GSM5745344", full.names = TRUE)
setwd("~/McGill/MethylVignettes/GDC Analysis/Controls")
for(i in Ctrls) getGEOSuppFiles(GEO = i, makeDirectory = T)
rm(i)
for(ct in Ctrls){
  f <- list.files(path = paste0("./Controls/", ct, collapse="/"), full.names = F)
  fldrName <- str_remove(f[1], "_Grn.idat.gz")
  file.rename(paste0("./Controls/", ct, collapse="/"), fldrName)
}

for(fs in f){
  proc <- list.files(fs, full.names = T)
  for(p in proc){ gunzip(p,remove = F) }
}

# Setting up actual datafiles.
# Extracting to top-level directories.
cases <- list.dirs("./cases")
cases <- dplyr::filter(as.data.frame(cases), grepl("maskedMethylation", cases))

# Sample sheet tsv -> csv.
tsvSS <- read_tsv("./cases/SingleFileCase/Sample_Sheet.tsv")
write.csv(tsvSS, file ="./cases/SingleFileCase/Sample_Sheet.csv")

# Bjarne Daenekas, Eilís Pérez, Fabio Boniolo, Sabina Stefan, Salvatore Benfatto, 
# Martin Sill, Dominik Sturm, David T W Jones, David Capper, Marc Zapatka, 
# Volker Hovestadt, Conumee 2.0: enhanced copy-number variation analysis from 
# DNA methylation arrays for humans and mice, Bioinformatics, Volume 40, Issue 2, 
# February 2024, btae029, https://doi.org/10.1093/bioinformatics/btae029

# Setting up controls.
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
x <- CNV.detail(x)
x <- CNV.segment(x)
segments <- CNV.write(x, what = "segments")
bins <- CNV.write(x, what = "bins")
write.csv(bins, "./Outputs/Conumee/PrzybylLabData/bins.csv")
for(segment in segments){
  write.csv(segment, paste0("./Outputs/Conumee/PrzybylLabData/", segment[["ID"]][1], ".csv"))
}

# Leiomyoma data.




# Mariani, M. P., Chen, J. A., Zhang, Z., Pike, S. C., & Salas, L. A. (2022). 
# MethylMasteR: A Comparison and Customization of Methylation-Based Copy Number 
# Variation Calling Software in Cancers Harboring Large Scale Chromosomal Deletions. 
# Frontiers in Bioinformatics, 2–2022. doi:10.3389/fbinf.2022.859828

# SeSAmE & MethylMasteR Custom.
input_dir <- "./cases/SingleFileCase"
output.dir <- paste0("./Outputs/MethylMaster/Sesame")
sample.sheet.path <- "./cases/SingleFileCase/Sample_Sheet.csv"
routine.run <- "sesame"

csvSS <- read.csv(sample.sheet.path)[,-c(1)]
write.csv(csvSS, sample.sheet.path)


methyl_master(
  routine                   = routine.run, #The routine to run
  input.dir                 = input_dir, #The input (idat.files) directory
  output.dir                = output.dir, #The output directory
  sample.sheet.path         = sample.sheet.path, #The path to the MethylMasteR sample sheet
  r.lib.path                = .libPaths()[1], #The path to the R Library path
  file.sep                  = "/", #For windows or "/" for Linux
  create.dir                = TRUE, #Whether to cretae directory if does not 
  #exist?
  save.seg                  = TRUE, #Whether to save segmentation results
  n.cores                   = 1, #Multicore does not work for all routines 
  #on all operating systems
  os.type                   = "linux", #Or "linux"
  proj                      = "TCGA-SARC", #"TCGA-BLCA" etc.
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
  sesame.data.cache         = "450K", #The default sesame reference platform
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
  olaps.split.field         = "File.ID", #Split field to ise during overlaps
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

# Tian Y, Morris TJ, Webster AP, Yang Z, Beck S, Andrew F, Teschendorff AE (2017). 
# “ChAMP: updated methylation analysis pipeline for Illumina BeadChips.” 
# Bioinformatics, btx513. doi:10.1093/bioinformatics/btx513. 
library(ChAMP)
data(EPICSimData)
CpG.GUI(arraytype="EPIC")
myNorm <- champ.norm(arraytype="EPIC")
myNorm_2 <- champ.norm(method="PBC",arraytype = "EPIC")
myNorm_3 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="SWAN",arraytype="EPIC")
myNorm_4 <- champ.norm(beta=myLoad_2$beta,mset=myLoad_2$mset,rgSet=myLoad_2$rgSet,method="FunctionalNormalization",arraytype="EPIC")
QC.GUI(arraytype="EPIC")
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
champ.SVD(beta=myNorm)
myDMP <- champ.DMP(arraytype = "EPIC")
DMP.GUI()
myDMR <- champ.DMR(arraytype = "EPIC")
myDMR_2 <- champ.DMR(arraytype = "EPIC",method="DMRcate",cores=1)
myDMR_3 <- champ.DMR(arraytype = "EPIC",method="ProbeLasso",compare.group=c("PrEC_cells","LNCaP_cells"))
DMR.GUI(DMR=myDMR,arraytype="EPIC",compare.group=c("PrEC_cells","LNCaP_cells"))
DMR.GUI(DMR=myDMR_2,arraytype="EPIC",compare.group=c("PrEC_cells","LNCaP_cells"))
DMR.GUI(DMR=myDMR_3,arraytype="EPIC",compare.group=c("PrEC_cells","LNCaP_cells"))
myGSEA <- champ.GSEA(DMP=myDMP[[2]],arraytype = "EPIC")
myBlock <- champ.Block(arraytype = "EPIC")
Block.GUI(arraytype="EPIC",compare.group=c("PrEC_cells","LNCaP_cells"))
myEpiMod <- champ.EpiMod(arraytype="EPIC")
myrefbase <- champ.refbase(arraytype = "EPIC")
myCNA <- champ.CNA(control = F,arraytype = "EPIC")

# Zhou, W., Triche, T. J., Jr, Laird, P. W., & Shen, H. (2018). 
# SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips
# in genomic deletions. Nucleic acids research, 46(20), e123. https://doi.org/10.1093/nar/gky691

library(sesame)
library(sesameData)
library(dplyr)
idat_dir <- "./LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM"
sdfs <- openSesame(idat_dir, func = NULL)
visualizeSegments(segs)




