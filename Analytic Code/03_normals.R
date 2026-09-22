### ============================================================
### 3. NORMALS
### ============================================================
# Baseline variability of each caller in normal-tissue controls.
# Writes: results/normals.csv

# Shared setup: libraries, helper functions, case lists / bin sizes.
# Run with the working directory set to the analysis root (the folder
# containing LabData/, Outputs/ and results/), since all data paths
# are built from getwd().
source(paste0(getwd(), "/Analytic Code/00_setup.R"))
source(paste0(getwd(), "/Analytic Code/01_functions.R"))
source(paste0(getwd(), "/Analytic Code/02_load_data.R"))

# Sanity-check normal-tissue controls: median/SD/min/max segment log2
# ratio per caller, across bin sizes, to confirm normals stay within
# an acceptable cutoff (i.e. show no spurious CNV calls).

normalsFrame10kb <- caseCorr(IDs = Normals, bin = 10000) %>% dplyr::mutate(bin = 10000)
normalsFrame100kb <- caseCorr(IDs = Normals, bin = 1e+05) %>% dplyr::mutate(bin = 1e+05)
normalsFrameDef <- caseCorr(IDs = Normals, bin = 50000) %>% dplyr::mutate(bin = 50000)
normalsFrame1Mb <- caseCorr(IDs = Normals, bin = 1e+06) %>% dplyr::mutate(bin = 1e+06)

NormalsFrame <- rbind(normalsFrame10kb, normalsFrame100kb, normalsFrameDef, 
                      normalsFrame1Mb) %>% arrange(bin) 
write.csv(NormalsFrame, file = paste0(getwd(), "/results/normals.csv"))

# We can observe none of the bins have anything outside of the cutoff

