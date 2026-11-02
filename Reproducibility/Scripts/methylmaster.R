# Imports
library(MethylMasteR)
library(sesameData)
sesameDataCacheAll()
library(tidyverse)

routine.run <- "custom"

# Both tumour cohorts are segmented against the five lab normal uterine smooth
# muscle samples (Sample_Group "normal"), via reference = "comparison" in the
# patched custom routine, instead of SeSAMe's internal EPIC.5.normal set.
#
# All IDATs are linked into one staging directory so both cohorts share a
# single input.dir (MethylMasteR points the ExperimentHub cache at input.dir).
# Each cohort's sample sheet selects its own tumours plus the same normals.
staging.dir <- file.path(tempdir(), "methylmaster_input")
dir.create(staging.dir, recursive = TRUE, showWarnings = FALSE)
idat.sources <- c(paste0(getwd(), "/LabData/MethylMasteRCombined"),  # LMS tumours + lab normals
                  paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files"))
for(src in idat.sources){
  for(f in list.files(src, pattern = "\\.idat$", full.names = TRUE)){
    dst <- file.path(staging.dir, basename(f))
    if(!file.exists(dst)) file.symlink(normalizePath(f), dst)
  }
}

lms.sheet <- read.csv(paste0(getwd(), "/LabData/MethylMasteRCombined/Sample_Sheet_Combined.csv"),
                      stringsAsFactors = FALSE)
lms.sheet <- lms.sheet[, names(lms.sheet) != "X"]
normals.sheet <- lms.sheet[lms.sheet$Sample_Group == "normal", ]

lm.ids <- read.csv(paste0(getwd(), "/LabData/LM_SNP_EPIC_array_data/EPIC_array_data_LM/idat_files/SampSheet.csv"),
                   stringsAsFactors = FALSE)
lm.names <- paste0(lm.ids$Sentrix_ID, "_", lm.ids$Sentrix_Position)
lm.sheet <- rbind(data.frame(Sample_Name = lm.names, Sample_Group = "tumor", Sample_Plate = NA,
                             Pool_ID = NA, Sentrix_Position = lm.ids$Sentrix_Position,
                             Sentrix_ID = lm.ids$Sentrix_ID, Basename = lm.names,
                             Platform = "EPIC", Batch = 1, stringsAsFactors = FALSE),
                  normals.sheet[, c("Sample_Name", "Sample_Group", "Sample_Plate", "Pool_ID",
                                    "Sentrix_Position", "Sentrix_ID", "Basename", "Platform", "Batch")])

sheet.dir <- file.path(tempdir(), "methylmaster_sheets")
dir.create(sheet.dir, showWarnings = FALSE)
write.csv(lms.sheet[, names(lm.sheet)], file.path(sheet.dir, "Sample_Sheet_LMS.csv"))
write.csv(lm.sheet, file.path(sheet.dir, "Sample_Sheet_LM.csv"))

cohorts <- list(
  LabLMS = list(sheet = file.path(sheet.dir, "Sample_Sheet_LMS.csv"),
                out   = paste0(getwd(), "/Outputs/MethylMaster/LabLMS")),
  LM     = list(sheet = file.path(sheet.dir, "Sample_Sheet_LM.csv"),
                out   = paste0(getwd(), "/Outputs/MethylMaster/LM"))
)

# Bin sizes (segmentation tile widths, bp) evaluated throughout the analysis.
# Output folders are named with R's default number formatting
# ("10000", "50000", "1e+05", "1e+06"), which is what the Analytic Code expects.
binLength <- c(10000, 50000, 1e+05, 1e+06)

mmstr <- function(input.dir, output.dir, sample.sheet.path, binSize){
  methyl_master(
    routine                   = routine.run, #The routine to run
    input.dir                 = input.dir, #The input (idat.files) directory
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
    proj                      = "TCGA-KIRC", #"TCGA-BLCA" etc.
    visualize                 = TRUE, #Whether to output plots,
    visualize.individual      = FALSE, #Whether to output indvidual sample plots,
    #only works for routine sesame
    reference                 = "comparison", #"comparison" or 'internal"; comparison uses the sample-sheet normals
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
    simplify.reduce           = weightedmean, #Equation to use during reduction
    custom.binsize            = binSize #Segmentation tile width (bp); exposed
    #by the patched sesame/MethylMasteR builds
  )
  
}

# Run MethylMasteR once per cohort and bin size, writing each run to its own
# folder (Outputs/MethylMaster/<cohort>/<bin>/<Sentrix>/autocorrected_regions.csv).
# output.dir is assigned globally because MethylMasteR's autocorrection step
# reads it from the calling environment.
# Optionally restrict to a subset of cohorts, e.g. MM_COHORTS=LM Rscript methylmaster.R
run.cohorts <- strsplit(Sys.getenv("MM_COHORTS", paste(names(cohorts), collapse = ",")), ",")[[1]]
stopifnot(all(run.cohorts %in% names(cohorts)))

for(cohort in run.cohorts){
  for(bin in binLength){
    output.dir <- paste0(cohorts[[cohort]]$out, "/", bin)
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    message("MethylMasteR: ", cohort, ", bin size ", bin, " -> ", output.dir)
    mmstr(input.dir = staging.dir,
          output.dir = output.dir,
          sample.sheet.path = cohorts[[cohort]]$sheet,
          binSize = bin)
  }
}
