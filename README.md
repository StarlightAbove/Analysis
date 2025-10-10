# Benchmarking detection of DNA copy number variation using high-density methylation microarrays
Codebase of the project developed over the summer at the Przybyl Lab, RI-MUHC. A key and guide to this repo and other details are being added to this README.


# Methods
This section is currently being written, and will cover the algorithms and visualizers being used.

**ALL CODEFILES ARE LOCATED IN ~/AnalysisCode**. It must be noted that there is a subdirectory within it known as "Old Analytical Code", which can be safely ignored as they were attempts to create a more natively functional algorithm without the use of GRanges.

Retrieval from the TCGA is not being covered, but that code is in ~/AnalysisCode/Old Analytical Code/integratedanalysis.R, and in the ~/NCI-GDC code folder, they are simple functions which do not require much explanation.

### Methylation to CNV Processing
All methylation to CNV processing occured within "~/AnalysisCode/labprocessing.R. There are three sections within that file, each looking at a different program. In this, there are three main variables to pay attention to:

***For Conumee2***
``` r
RGset <- read.metharray.exp("./LabData/LMS_SNP_EPIC_array_data", recursive = T, verbose = T)
MSet <- preprocessIllumina(RGset)
CSet <- read.metharray.exp("./LabData/Normal_smooth_muscle_EPIC_data", recursive = T, verbose = T)
CSet <- preprocessIllumina(CSet)
anno <- CNV.create_anno(array_type = "EPIC")
load.data <- CNV.load(MSet)
load.controls <- CNV.load(CSet)
```
The sets of control and lab data were loaded, preprocessed, and an annotation was created. Eventually, these were loaded to the typical locations for data and its controls.

***For SeSAMe***
``` r
idat_dir <- "./LabData/LMS_SNP_EPIC_array_data/EPIC_array_data_LMS/idat_files"
sdfs <- openSesame(idat_dir, func = NULL)
```
SeSAMe internally creates its own idat references, so only a reference to the directory directly was given, and then opened.

***For MethylMasteR***
``` r
...
input.dir <- paste0(getwd(), "/LabData/MethylMasteRCombined")
output.dir <- paste0(getwd(),"/Outputs/MethylMaster/LabLMS")
routine.run <- "custom"
sample.sheet.path <- paste0(getwd(), "/LabData/MethylMasteRCombined/Sample_Sheet_Combined.csv")
...
```
input.dir is the input files of the .idats for the program. output.dir is where the outputs should go, routine.run is the analytical routine, and sample.sheet.path reads the sample sheet of the .idats, which were custom created to match the requirements of the program. It must be noted that lab controls were used, just in the main body of the command.
The rest of the code came directly from their vignettes, with some basic for-loops which were trivially used to organize output. These files can ALL be found in ~/Outputs, organized by program.

### Preprocessing for Accuracy
The prep.rocessing system is rather simple, and spread across three different codefiles in the AnalysisCode folder, specifically grangesSesame.R, grangesConumee.R and grangesMethyl.R, each refering to their respective software.
Preprocessing files take output from the code as well as its SNP files, and then places them into a single file listing everything, while labelling each in a singular internal reference. Conumee's way of referring to its outputs was taken as the standard.
The main problem we would run into here was that each file, the LM, LMS and laboratory LMS files had different "unique" identifiers, so therefore, they demanded different retrieval techniques. These can be seen in the labelling in the tables below. 

***Lab LMS***:

This is the one piece of code where all three different programs were correlated by STT, therefore there is only one function. The comments added to the code below explain it.
``` r
labLMSProc <- function(STTq, Technology){
  # Correlate by STT information between methylation Sentrix and SNP data.

  # The correlation sheet here is a csv that was generated which lists every file with its STT, therefore by filtering for STT, we can get the files it correlates to.
  correlationSheet <- read.csv("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/correlative.csv") %>% filter(STT == STTq)

  # From the correlation sheet, it can take the IDs for the CNV as well as the methylation data via Sentrix Location and Position, and retrieve all of the files automatically.
  cnvMatch <- read.csv(paste0("~/Work/Analysis/LabData/LMS_SNP_EPIC_array_data/SNP_array_data_LMS/CNV_calls/CNVCallsCSV/", 
                              correlationSheet$CNV_Label[1], "_events.csv"))
  
  methylMatch <- read.csv(paste0("~/Work/Analysis/Outputs/MethylMaster/LabLMS/", 
                                 paste0(correlationSheet$Sentrix_ID[1],"_", 
                                        correlationSheet$Sentrix_Position[1], "/"), 
                                 "autocorrected_regions.csv"))
  
  conumeeMatch <- read.csv(paste0("~/Work/Analysis/Outputs/Conumee/LabLMS/", paste0(correlationSheet$Sentrix_ID[1],"_", 
                                                                                    correlationSheet$Sentrix_Position[1],".csv")))
  
  sesameMatch <- read.csv(paste0("~/Work/Analysis/Outputs/SeSAMe/LabLMS/", paste0("segments_",correlationSheet$Sentrix_ID[1],"_", 
                                                                                   correlationSheet$Sentrix_Position[1],".csv")))

  # The CNV files are renamed to the internal reference, with loc.start indicating the start, loc.end indicating the end of the region, seg.mean indicating the mean,
  # chrom indicating the specific chromosome, and type indicating which software it is from.

  cnvMatch <- tidyr::separate_wider_delim(tidyr::separate_wider_delim(cnvMatch, `Chromosome.Region`, names = c("chrom", "loc"), delim = ":"), `loc`,
                                          names = c("loc.start", "loc.end"), delim = "-") %>% mutate(type = "SNP") %>% dplyr::rename(CNVStatus = "Event") %>%
    mutate(CNVStatus = case_when(
      CNVStatus == "CN Loss" ~ "Deletion", 
      CNVStatus == "CN Gain" ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% dplyr::rename(seg.mean = "Probe.Median") %>% 
    dplyr::select(c(chrom, loc.start, loc.end, CNVStatus, seg.mean, type)) %>% mutate(loc.start = as.numeric(gsub(",","",loc.start)), loc.end = as.numeric(gsub(",","",loc.end)))
  cnvMatch$chrom <- as.numeric(str_replace_all(cnvMatch$chrom, "chr", ""))
  cnvMatch <- cnvMatch %>% filter(!is.na(chrom))

  # Each technology gives its outputs with different columns, which are then renamed based on the technology, and then added to the "combinedset" to be outputted.
  if(Technology == "MethylMaster") {
    methylMatch <- methylMatch %>% dplyr::rename(
      seg.mean = "Mean", loc.start = "bp.Start", 
      loc.end = "bp.End") %>% dplyr::select(
        seg.mean, loc.start, loc.end, Chromosome) %>% mutate(
          type = "MethylMaster") %>% dplyr::rename(chrom = "Chromosome") %>% 
      dplyr::mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% filter(!is.na(chrom)) %>%
      dplyr::mutate(loc.start = as.numeric(loc.start), loc.end = as.numeric(loc.end), seg.mean = as.numeric(seg.mean)) %>% mutate(CNVStatus = case_when(
        seg.mean <= -0.2 ~ "Deletion", 
        seg.mean >= 0.2 ~ "Amplification", 
        TRUE ~ "Normal"
      ))
    combinedSet <- rbind(methylMatch, cnvMatch) %>% arrange(chrom)
  }
  
  if(Technology == "Conumee"){
    conumeeMatch <- conumeeMatch %>% dplyr::select(chrom, loc.start, loc.end, seg.mean) %>% mutate(CNVStatus = case_when(
      seg.mean <= -0.2 ~ "Deletion", 
      seg.mean >= 0.2 ~ "Amplification", 
      TRUE ~ "Normal"
    )) %>% mutate(chrom = as.numeric(gsub("chr", "", chrom))) %>% mutate(type = "Conumee")
    combinedSet <- rbind(conumeeMatch, cnvMatch) %>% arrange(chrom)
  }
  
  if(Technology == "Sesame"){
    sesameMatch <- sesameMatch %>% dplyr::select(c("chrom", "loc.start", 
                                                   "loc.end", "seg.mean")) %>% dplyr::mutate(
                                                     chrom = str_remove_all(chrom, "chr")) %>% filter(
                                                       !(chrom == "X") & !(chrom == "Y")) %>% mutate(
                                                         chrom = as.numeric(chrom)) %>% arrange(chrom) %>% mutate(
                                                           CNVStatus = case_when(seg.mean > 0.3 ~ "Amplification",
                                                                                 seg.mean < -0.3 ~ "Deletion",
                                                                                 TRUE ~ "Normal"), type = "SeSAMe") 
    combinedSet <- rbind(sesameMatch, cnvMatch) %>% arrange(chrom)
  }
  
  
  
  combinedSet
}
```



### Accuracy Determination

### Pearson Correlation Coefficient Computation

### Visualization
NOTE: Due to recent Visualization changes, there is a slight artifacting in the SNP Array line within Manhattan plots. These will be changed and swapped as soon as possible. This was due to the addition of "Mixedsort" to a labelling system, which was not later removed.

# Indexes for all statistics
Details the internal sample references, their correlated identifiers, and the filepath of the Manhattan plot.

### "LMS" - LMS cases from TCGA-provided samples.
| Sample Number (Internal) | TCGA Identifier | File Path                                    | Most/Least Accurate? | Comments |
|--------------------------|-----------------|----------------------------------------------|----------------------|----------|
| Sample 1                 | TCGA-3B-A9HQ    | ./Statistics/LMS/AllSamples/TCGA-3B-A9HQ.png |                      |          |
| Sample 2                 | TCGA-3B-A9HX    | ./Statistics/LMS/AllSamples/TCGA-3B-A9HX.png |                      |          |
| Sample 3                 | TCGA-DX-A3UC    | ./Statistics/LMS/AllSamples/TCGA-DX-A3UC.png |                      |          |
| Sample 4                 | TCGA-DX-A48O    | ./Statistics/LMS/AllSamples/TCGA-DX-A48O.png |                      |          |
| Sample 5                 | TCGA-DX-A48P    | ./Statistics/LMS/AllSamples/TCGA-DX-A48P.png |                      |          |
| Sample 6                 | TCGA-DX-A6BA    | ./Statistics/LMS/AllSamples/TCGA-DX-A6BA.png | ✖ Least Accurate     | Labelled as lAcc.png. |
| Sample 7                 | TCGA-DX-A6Z2    | ./Statistics/LMS/AllSamples/TCGA-DX-A6Z2.png |                      |          |
| Sample 8                 | TCGA-DX-A7EL    | ./Statistics/LMS/AllSamples/TCGA-DX-A7EL.png |                      |          |
| Sample 9                 | TCGA-DX-A7EN    | ./Statistics/LMS/AllSamples/TCGA-DX-A7EN.png | ✔ Most Accurate      | Labelled as mAcc.png  |
| Sample 10                | TCGA-K1-A42W    | ./Statistics/LMS/AllSamples/TCGA-K1-A42W.png |                      |          |

### "LabLMS" - LMS cases from laboratory-provided samples.
| Sample Number (Internal) | STT Identifier | Filepath to Manhattan Plot              | Most/Least Accurate | Comments |
|--------------------------|----------------|-----------------------------------------|---------------------|----------|
| Sample 1                 | 9202           | ./Statistics/LabLMS/AllSamples/9202.png |                     |          |
| Sample 2                 | 9203           | ./Statistics/LabLMS/AllSamples/9203.png | ✔ Most Accurate     |          |
| Sample 3                 | 9327           | ./Statistics/LabLMS/AllSamples/9327.png |                     |          |
| Sample 4                 | 9328           | ./Statistics/LabLMS/AllSamples/9328.png |                     |          |
| Sample 5                 | 9337           | ./Statistics/LabLMS/AllSamples/9337.png |                     |          |
| Sample 6                 | 9338           | ./Statistics/LabLMS/AllSamples/9338.png |                     |          |
| Sample 7                 | 9350           | ./Statistics/LabLMS/AllSamples/9350.png |                     |          |
| Sample 8                 | 9353           | ./Statistics/LabLMS/AllSamples/9353.png | ✖ Least Accurate    |          |
| Sample 9                 | 9354           | ./Statistics/LabLMS/AllSamples/9354.png |                     |          |
| Sample 10                | 9355           | ./Statistics/LabLMS/AllSamples/9355.png |                     |          |
| Sample 11                | 9356           | ./Statistics/LabLMS/AllSamples/9356.png |                     |          |
| Sample 12                | 9357           | ./Statistics/LabLMS/AllSamples/9357.png |                     |          |
| Sample 13                | 9358           | ./Statistics/LabLMS/AllSamples/9358.png |                     |          |

### "LM" - LM cases from laboratory-provided samples.
| Sample Number (Internal) | Sentrix Identifier  | Filepath to Manhattan Plot                         | Most/Least Accurate | Comments |
|--------------------------|---------------------|----------------------------------------------------|---------------------|----------|
| Sample 1                 | 203219650162_R08C01 | ./Statistics/LM/AllSamples/203219650162_R08C01.png |                     |          |
| Sample 2                 | 203219650188_R08C01 | ./Statistics/LM/AllSamples/203219650188_R08C01.png | ✖ Least Accurate    |          |
| Sample 3                 | 203219750019_R05C01 | ./Statistics/LM/AllSamples/203219750019_R05C01.png |                     |          |
| Sample 4                 | 203219750053_R05C01 | ./Statistics/LM/AllSamples/203219750053_R05C01.png |                     |          |
| Sample 5                 | 203219760009_R08C01 | ./Statistics/LM/AllSamples/203219760009_R08C01.png |                     |          |
| Sample 6                 | 203220070150_R02C01 | ./Statistics/LM/AllSamples/203220070150_R02C01.png |                     |          |
| Sample 7                 | 203225140266_R03C01 | ./Statistics/LM/AllSamples/203225140266_R03C01.png |                     |          |
| Sample 8                 | 203257030015_R06C01 | ./Statistics/LM/AllSamples/203257030015_R06C01.png |                     |          |
| Sample 9                 | 203257030060_R08C01 | ./Statistics/LM/AllSamples/203257030060_R08C01.png | ✔ Most Accurate     |          |

Most accurate Manhattan plots are usually labelled as mAcc.png, and least accurate as lAcc.png.

