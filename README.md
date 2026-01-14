# Benchmarking detection of DNA copy number variation using high-density methylation microarrays

Codebase of the project developed over the summer at the Przybyl Lab, RI-MUHC. A key and guide to this repo and other details are being added to this README.

## News

-   **New Commit** (29/12/2025, 19:00 IST): Added sections towards the end of this document in relation to new heatmap graphs, genomic labeling, as well as CEL processing.
    Check commits bb4b265c8f5c2e99785f50b4278cb60ed6648bb1 and forward till 01/01/2026 for more details and tracking.

# Methods

This section is currently being written, and will cover the algorithms and visualizers being used.

**ALL CODEFILES ARE LOCATED IN \~/AnalysisCode**.
It must be noted that there is a subdirectory within it known as "Old Analytical Code", which can be safely ignored as they were attempts to create a more natively functional algorithm without the use of GRanges.

Retrieval from the TCGA is not being covered, but that code is in \~/AnalysisCode/Old Analytical Code/integratedanalysis.R, and in the \~/NCI-GDC code folder, they are simple functions which do not require much explanation.

### Methylation to CNV Processing

All methylation to CNV processing occured within "\~/AnalysisCode/labprocessing.R. There are three sections within that file, each looking at a different program. In this, there are three main variables to pay attention to:

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

The sets of control and lab data were loaded, preprocessed, and an annotation was created.
Eventually, these were loaded to the typical locations for data and its controls.

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

input.dir is the input files of the .idats for the program.
output.dir is where the outputs should go, routine.run is the analytical routine, and sample.sheet.path reads the sample sheet of the .idats, which were custom created to match the requirements of the program.
It must be noted that lab controls were used, just in the main body of the command.
The rest of the code came directly from their vignettes, with some basic for-loops which were trivially used to organize output.
These files can ALL be found in \~/Outputs, organized by program.

### Preprocessing for Accuracy

The prep.rocessing system is rather simple, and spread across three different codefiles in the AnalysisCode folder, specifically grangesSesame.R, grangesConumee.R and grangesMethyl.R, each refering to their respective software.
Preprocessing files take output from the code as well as its SNP files, and then places them into a single file listing everything, while labelling each in a singular internal reference.
Conumee's way of referring to its outputs was taken as the standard.
The main problem we would run into here was that each file, the LM, LMS and laboratory LMS files had different "unique" identifiers, so therefore, they demanded different retrieval techniques.
These can be seen in the labelling in the tables below.

***Lab LMS***:

This is the one piece of code where all three different programs were correlated by STT, therefore there is only one function.
The comments added to the code below explain it.

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

***LMS and LM files***: LMS and LM files are split across two functions, just due to the idea fo consolidation at that moment not popping into my head.
They are very similar however, and it is a future goal in the code cleaning to get that done.
In short, they work the same way as the code above, which is why I am choosing not to go into detail.
They either use the TCGA identifier, or the Sentrix identifier directly to retrieve data, and rename it into the internal reference which I mentioned above, returning a combined data.frame, with both CNV calls from SNP arrays and from software, all labelled with their "type".

An output from this code is sampled below.

``` r
chrom loc.start   loc.end seg.mean CNVStatus     type   
    <chr>     <dbl>     <dbl>    <dbl> <chr>         <chr>  
  1 1        599171   3185375  -0.0629 Normal        SNP    
  2 1       3186748   3906704  -0.753  Deletion      SNP    
  3 1       3974037  11900776   0.033  Normal        SNP    
  4 1      11909975  16482199   0.492  Amplification SNP    
  5 1      16487426  17040250   0.261  Normal        SNP    
  6 1      17040376  25256850   0.488  Amplification SNP    
  7 1      25336819  31663315   0.485  Amplification SNP    
  8 1      31670172  31972468   0.022  Normal        SNP    
  9 1      31972914  45785677   0.491  Amplification SNP    
 10 1      45903822  48342115   0.0456 Normal        SNP    
 11 1      48343508  50572289   0.485  Amplification SNP    
 12 1      50574582  53957577   0.0636 Normal        SNP    
 13 1      54137895  58254237   0.0481 Normal        SNP
# There of course exist software outputs, but these are sorted across by type.
```

### Accuracy Determination

For accuracy determination, a dataframe of CNV-methylation correlated is required.
This then goes through either "fpCheck.R" or "gRangeAccModel.R".
fpCheck is specifically used for LM data due to how narrow the changes are.
Both of them have generally the same code, with fpCheck having the input of the HG19 reference genome to compute across the whole genome.

***Analysing gRangeAccModel.R***: The entire model is located in a function named "AccuracyModel".

``` r
pred_df <- t1 %>% filter(!(type == "SNP"))
  truth_df <- t1 %>% filter(type == "SNP")
  
  truth_gr <- GRanges(seqnames = truth_df$chrom,
                      ranges = IRanges(start = truth_df$loc.start, end = truth_df$loc.end),
                      CNV = truth_df$CNVStatus)
  
  pred_gr <- GRanges(seqnames = pred_df$chrom,
                     ranges = IRanges(start = pred_df$loc.start, end = pred_df$loc.end),
                     CNV = pred_df$CNVStatus)
```

This is simple, it takes the table given to it, and sorts it into software CNV calls and CNV calls from SNP arrays, and converts them into gRanges.
IRanges are used for the location of the CNV since it will simplify the math later on.

``` r
hits <- findOverlaps(pred_gr, truth_gr)
overlap_ranges <- pintersect(pred_gr[queryHits(hits)], truth_gr[subjectHits(hits)])
```

These two commands are the "meat" of the algorithm.
"hits" stores the findOverlaps output, which is specifically an IRanges function.
This function will look like this:

``` r
Hits object with 248 hits and 0 metadata columns:
        queryHits subjectHits
        <integer>   <integer>
    [1]         1           1
    [2]         2           1
    [3]         2           2
    [4]         3           3
    [5]         3           4
    ...       ...         ...
  [244]        60         212
  [245]        60         213
  [246]        60         214
  [247]        60         215
  [248]        60         216
```

What this essentially is doing is that it is numbering each instance in the software output, and then matching it based on its range to the "subject", in this case, the CNVs called from the SNPs.
This works both ways, so there is no need to consider any reverse or forward implication.

pIntersect then takes the results from this, to "intersect" parallel, essentially, imagine that each of these were groups by their numbers, and then the intersect function from set theory was applied.
That, in a nutshell, is pintersect.
This output will look like this:

``` r
GRanges object with 248 ranges and 2 metadata columns:
        seqnames            ranges strand |         CNV       hit
           <Rle>         <IRanges>  <Rle> | <character> <logical>
    [1]        1    599171-3075000      * |      Normal      TRUE
    [2]        1   3125000-3185375      * |    Deletion      TRUE
    [3]        1   3186748-3822634      * |    Deletion      TRUE
    [4]        1  4072634-11900776      * |      Normal      TRUE
    [5]        1 11909975-16482199      * |      Normal      TRUE
    ...      ...               ...    ... .         ...       ...
  [244]       22 16373925-18822597      * |      Normal      TRUE
  [245]       22 18917479-23432731      * |      Normal      TRUE
  [246]       22 23433022-42496590      * |      Normal      TRUE
  [247]       22 42546956-49004270      * |      Normal      TRUE
  [248]       22 49004572-50796027      * |      Normal      TRUE
```

It shows the different IRanges along which these intersect, as well as their CNV character.
seqnames, of course, correlates to chromosome.
Now, all that is left to compute their sizes and average them out to complete the dataset (that's the hard part done!).

``` r
overlap_df <- data.frame(
    Chromosome = as.character(seqnames(overlap_ranges)),
    width = width(overlap_ranges),
    pred_cnv = mcols(pred_gr)$CNV[queryHits(hits)],
    truth_cnv = mcols(truth_gr)$CNV[subjectHits(hits)]
  )
  
  weighted_matrix <- overlap_df %>%
    group_by(truth_cnv, pred_cnv) %>%
    summarise(total_bp = sum(width), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pred_cnv, values_from = total_bp, values_fill = 0)
```

overlap_df is merely consolidating all the different data points into the values which we need.
Chromosome is of course seqnames, width returns the difference of each IRange, and pred_cnv is taking the "mcols" (the values beyond the \| line, and then looking for their CNV which correlates to their "hits", and the same is done from the truth_cnv. The output is given below:

``` r
    Chromosome    width pred_cnv     truth_cnv
1            1  2475830   Normal        Normal
2            1    60376 Deletion        Normal
3            1   635887 Deletion      Deletion
4            1  7828143   Normal        Normal
5            1  4572225   Normal Amplification
6            1   552825   Normal        Normal
7            1  8216475   Normal Amplification
8            1  6326497   Normal Amplification
9            1   302297   Normal        Normal
10           1 13812764   Normal Amplification
...
```

Now, we have to find the matrix of different values, think of this as a confusion matrix.
This is what the weighted_matrix command does, which makes it a lot simpler for use across the entire sample.
It groups by the different predictions, sums their "widths", and then drops the grouping.
It then pivots it wider to show the confusion matrix.
The output is below:

``` r
  truth_cnv      Deletion     Normal
  <chr>             <int>      <int>
1 Amplification  11622220  308199672
2 Deletion      381256807   50694212
3 Normal        308760825 1687770847
```

Note that this matrix has been bypassed in the latest iteration, but can give interesting insight, so I have left it in the code.
This can then be used to compute the total "truth" vs the total basepairs covered in the first place.

``` r
TP_weighted <- sum(overlap_df$width[overlap_df$truth_cnv == overlap_df$pred_cnv])
Total_weighted <- sum(overlap_df$width)
weighted_accuracy <- TP_weighted / Total_weighted
```

This is done rather simply, TP_weighted takes the widths which match, and the Total_weighted sums the whole thing.
Then the ratio of the two is taken.
Per chromosome, the same principle is used, but all that is done is grouped by chromosome before computation:

``` r
accuracy_by_chr <- overlap_df %>%
    group_by(Chromosome) %>%
    summarise(
      TP_bp = sum(width[pred_cnv == truth_cnv]),
      Total_bp = sum(width),
      Accuracy = TP_bp / Total_bp
    )
```

This all leads to the final output as given below:

``` r
[[1]]
# A tibble: 22 × 4
   Chromosome     TP_bp  Total_bp Accuracy
   <chr>          <int>     <int>    <dbl>
 1 1          197517883 246108941    0.803
 2 10         133185948 133185948    1    
 3 11         103276264 132148541    0.782
 4 12          96383261 129706240    0.743
 5 13          94543856  94543856    1    
 6 14          76639185  87229477    0.879
 7 15          19122580  79086734    0.242
 8 16          89674915  89972386    0.997
 9 17          40235548  80076561    0.502
10 18          77886221  77886221    1    
# ℹ 12 more rows
# ℹ Use `print(n = ...)` to see more rows

[[2]]
[1] 0.7528378
```

fpCheck works exactly on these same lines with some important corrections.
Up until the overlaps, it is the exact same, and then it continues on an extended subroutine.

``` r
  non_overlapping_pred <- pred_gr[-queryHits(hits)]
  fp_df <- as.data.frame(non_overlapping_pred) %>%
    filter(CNV != "Normal") %>%
    group_by(seqnames) %>%
    summarise(FP_bp = sum(width), .groups = "drop") %>% dplyr::rename(Chromosome = "seqnames")
  print(fp_df)
```

It computes what is not within the queryHits (the truth) but in pred_gr (the software CNV).
It then computes a data.frame where it filters out the normals, and then groups by the chromosomes, and computes the false positive basepairs finally.

``` r
  non_overlapping_truth <- truth_gr[-subjectHits(hits)]
  fn_df <- as.data.frame(non_overlapping_truth) %>%
    filter(CNV != "Normal") %>%
    group_by(seqnames) %>%
    summarise(FN_bp = sum(width), .groups = "drop") %>% dplyr::rename(Chromosome = "seqnames")
  print(fn_df)
```

It computes what is not within the subjectHits (the software CNV) but in truth_gr (the truth set).
It then computes a data.frame where it filters out the normals, and then groups by the chromosomes, and computes the false negative basepairs finally.

``` r
  final_eval <- hg19_chr_sizes %>%
    left_join(tp_df, by = "Chromosome") %>%
    left_join(fp_df, by = "Chromosome") %>%
    left_join(fn_df, by = "Chromosome") %>%
    mutate(across(c(TP_bp, FP_bp, FN_bp), ~replace_na(., 0))) %>%
    mutate(
      CNV_bp = TP_bp + FP_bp + FN_bp,
      TN_bp = Genome_bp - CNV_bp,
      Accuracy = (TP_bp + TN_bp) / Genome_bp,
      CNV_Only_Accuracy = ifelse((TP_bp + FP_bp + FN_bp) == 0, NA,
                                 TP_bp / (TP_bp + FP_bp + FN_bp))
    ) %>%
    arrange(as.numeric(Chromosome))
  
  accuracy <- sum(final_eval$Accuracy) / 22
```

All of this is then merged together, with all the NA replaced by a mutate, then the CNV basepairs are computed, with true negative and positive giving the accuracy across the entire genomic reference of HG19.
The rest is within the CNVs only, such as the method in gRangeAccModel.
This is of course, computed across each chromosome, so the sum of the final accuracies divided by 22 gives the final sample-wide accuracy.

This gives us some extra information and seemingly gave more sensible results on LM.
However, this is a lot more unreliable and gave some \>1 accuracies, therefore I do not trust it as much.

### Pearson Correlation Coefficient Computation (pcc.R)

This hijacks a part of the normal accuracy algorithm to compute a Pearson Correlation Coefficient.
Similar the above code, until overlaps_df, the code is the same.
The only difference is the further computation of accuracy is removed and instead replaced by:

``` r
  cor(overlap_df$pred_cnv, overlap_df$truth_cnv, method = "pearson")
```

This will compute the statuses over each other, and acts as a state variable with three different modes, whcih then computes for their accurate correlation.

### Visualization (visualize.R and consolidation.R)

NOTE: Due to recent Visualization changes, there is a slight artifacting in the SNP Array line within Manhattan plots for Lab LMS. I am currently analysing them and will try to get them fixed ASAP.

consolidation.R contains all the files which synthesize statistics as well as the Manhattan plots.
The Manhattan genomic plots are derived from visualize.R.
These use the outputs from post-preprocessing, with all software integrated into one file, and that is merely passed through to ggplot2.
The implementation is rather simplistic and will be described below:

``` r
# Calculate chromosome cumulative positions
  chr_lengths <- df %>%
    group_by(chrom) %>%
    summarize(chr_len = max(loc.end), .groups = "drop") %>%
    arrange(chrom) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
    mutate(chr_mid = chr_start + chr_len / 2)
  # Join to get cumulative start and end positions
  df <- df %>%
    left_join(chr_lengths, by = "chrom") %>%
    mutate(
      start_cum = loc.start + chr_start,
      end_cum = loc.end + chr_start
    ) %>%mutate(type = ifelse(type == "SNP", "SNP Array", type)) %>% mutate(type = as.factor(type))
```

It first calculates the chromosomal lengths, and then finds the middle.
It then adds the actual start and end linearly to the actual dataframe that is then visualized.

``` r
chr_boundaries <- chr_lengths %>%
    mutate(x = chr_start) %>%
    select(chrom, x)
  
  x_breaks <- chr_lengths$chr_mid
  x_labels <- paste0("chr", chr_lengths$chrom)
```

The chromosome limits are found, and then the chromosomes and x (the chromosome start are selected).
The lengths are also added as "breaks", giving the vertical bars.
The x-labels allow us to mention which chromosome is being represented on the linear graph.

``` r
# Plot
  p <- ggplot(df, aes(x = start_cum, xend = end_cum, y = seg.mean, yend = seg.mean)) +
    geom_segment(aes(color = type), size = 0.7, alpha = 0.8) +
    geom_vline(data = chr_boundaries, aes(xintercept = x), color = "grey70", linetype = "dashed") +
    # scale_color_manual(values = c("Amplification" = "red", "Deletion" = "blue", "Normal" = "black")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("SeSAMe" = "red", "MethylMaster" = "blue", "Conumee" = "green", "SNP Array" = "black")) +
    labs(
      x = "Genomic Position (across chromosomes)",
      y = "Segment Mean (log2 ratio)",
      title = "CNV Segments Across Genome"
    ) + geom_hline(yintercept = -0.2, linetype = "dotted", color = "black") + 
    geom_hline(yintercept = 0.2, linetype = "dotted", color = "black") + 
    theme_minimal() +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.major.x = element_blank(),
      legend.position = "bottom"
    )
```

It is then plotted, with the x being from the cumulative start, to the cumulative end, and the same for the y.
Then, segments are created by type, and vertical lines are added to segment chromosomes.
Then, the breaks for the vertical lines are added, and a scale color is added so that each software and SNP array can be represented.
Then, x and y labels at large are added, and horizontal lines are added to represent the limits used for the deletion-amplification inferrence.
Finally, a theme is set, and "p" is returned to produce the plot.

# Indexes for all statistics

Details the internal sample references, their correlated identifiers, and the filepath of the Manhattan plot.

### "LMS" - LMS cases from TCGA-provided samples.

| Sample Number (Internal) | TCGA Identifier | File Path | Most/Least Accurate? | Comments |
|-------------|-------------|----------------------|-------------|-------------|
| Sample 1 | TCGA-3B-A9HQ | ./Statistics/LMS/AllSamples/TCGA-3B-A9HQ.png |  |  |
| Sample 2 | TCGA-3B-A9HX | ./Statistics/LMS/AllSamples/TCGA-3B-A9HX.png |  |  |
| Sample 3 | TCGA-DX-A3UC | ./Statistics/LMS/AllSamples/TCGA-DX-A3UC.png |  |  |
| Sample 4 | TCGA-DX-A48O | ./Statistics/LMS/AllSamples/TCGA-DX-A48O.png |  |  |
| Sample 5 | TCGA-DX-A48P | ./Statistics/LMS/AllSamples/TCGA-DX-A48P.png |  |  |
| Sample 6 | TCGA-DX-A6BA | ./Statistics/LMS/AllSamples/TCGA-DX-A6BA.png | ✖ Least Accurate | Labelled as lAcc.png. |
| Sample 7 | TCGA-DX-A6Z2 | ./Statistics/LMS/AllSamples/TCGA-DX-A6Z2.png |  |  |
| Sample 8 | TCGA-DX-A7EL | ./Statistics/LMS/AllSamples/TCGA-DX-A7EL.png |  |  |
| Sample 9 | TCGA-DX-A7EN | ./Statistics/LMS/AllSamples/TCGA-DX-A7EN.png | ✔ Most Accurate | Labelled as mAcc.png |
| Sample 10 | TCGA-K1-A42W | ./Statistics/LMS/AllSamples/TCGA-K1-A42W.png |  |  |

### "LabLMS" - LMS cases from laboratory-provided samples.

| Sample Number (Internal) | STT Identifier | Filepath to Manhattan Plot | Most/Least Accurate | Comments |
|-------------|-------------|--------------------|-------------|-------------|
| Sample 1 | 9202 | ./Statistics/LabLMS/AllSamples/9202.png |  |  |
| Sample 2 | 9203 | ./Statistics/LabLMS/AllSamples/9203.png | ✔ Most Accurate |  |
| Sample 3 | 9327 | ./Statistics/LabLMS/AllSamples/9327.png |  |  |
| Sample 4 | 9328 | ./Statistics/LabLMS/AllSamples/9328.png |  |  |
| Sample 5 | 9337 | ./Statistics/LabLMS/AllSamples/9337.png |  |  |
| Sample 6 | 9338 | ./Statistics/LabLMS/AllSamples/9338.png |  |  |
| Sample 7 | 9350 | ./Statistics/LabLMS/AllSamples/9350.png |  |  |
| Sample 8 | 9353 | ./Statistics/LabLMS/AllSamples/9353.png | ✖ Least Accurate |  |
| Sample 9 | 9354 | ./Statistics/LabLMS/AllSamples/9354.png |  |  |
| Sample 10 | 9355 | ./Statistics/LabLMS/AllSamples/9355.png |  |  |
| Sample 11 | 9356 | ./Statistics/LabLMS/AllSamples/9356.png |  |  |
| Sample 12 | 9357 | ./Statistics/LabLMS/AllSamples/9357.png |  |  |
| Sample 13 | 9358 | ./Statistics/LabLMS/AllSamples/9358.png |  |  |

### "LM" - LM cases from laboratory-provided samples.

| Sample Number (Internal) | Sentrix Identifier | Filepath to Manhattan Plot | Most/Least Accurate | Comments |
|-------------|-------------|-----------------------|-------------|-------------|
| Sample 1 | 203219650162_R08C01 | ./Statistics/LM/AllSamples/203219650162_R08C01.png |  |  |
| Sample 2 | 203219650188_R08C01 | ./Statistics/LM/AllSamples/203219650188_R08C01.png | ✖ Least Accurate |  |
| Sample 3 | 203219750019_R05C01 | ./Statistics/LM/AllSamples/203219750019_R05C01.png |  |  |
| Sample 4 | 203219750053_R05C01 | ./Statistics/LM/AllSamples/203219750053_R05C01.png |  |  |
| Sample 5 | 203219760009_R08C01 | ./Statistics/LM/AllSamples/203219760009_R08C01.png |  |  |
| Sample 6 | 203220070150_R02C01 | ./Statistics/LM/AllSamples/203220070150_R02C01.png |  |  |
| Sample 7 | 203225140266_R03C01 | ./Statistics/LM/AllSamples/203225140266_R03C01.png |  |  |
| Sample 8 | 203257030015_R06C01 | ./Statistics/LM/AllSamples/203257030015_R06C01.png |  |  |
| Sample 9 | 203257030060_R08C01 | ./Statistics/LM/AllSamples/203257030060_R08C01.png | ✔ Most Accurate |  |

Most accurate Manhattan plots are usually labelled as mAcc.png, and least accurate as lAcc.png.

# Binning Systems in Conumee 2.0, SeSAMe & MethylMasteR.

## LabLMS

Graphs Path: Analysis/Statistics/LabLMS/byBin It is organized by the binsize, and then also the case that it comes from.

## LMS

Graphs path: currently under computation.
Holding off till MethylMasteR stable.

## LM

Graphs path: currently under computation.
Holding off till MethylMasteR stable.

## Normals

Normals have been computed and exists in: - Analysis/Outputs/Conumee/LabNormals: organized by binsize.
- Analysis/Outputs/SeSAMe/Normals: organized by binsize.
The default ones are in 50000.
The graphs have been computed and can be found at: Analysis/Statistics/Normals/byBin, organized by bin size, with the default being 50000.

# CEL Processing

Located at: /Analysis/LabData/ReprocessedLMS for LMS data (Laboratory) and /Analysis/LabData/ReprocessedLM for LM data.
Currently, due to new information, B-allele data is being integrated into the analysis in order to further refine data-collection from ASCAT & SEQUENZA, and to solve issues related to over-calling of certain genotypes. This approach is based upon Connolly et al., 2014 & Wang, 2007. 

Things to note: The **correlation_w_FASST2.csv** file is the core characterization of this new data and shows its alignment with FASST2.
If you are trying to evaluate alignment, please look for that data first.
It characterizes both the data from ASCAT as well as the SEQUENZA methods as offered by EaCoN.

# Genomic Plots + Heatmap

The genomic plots have now been combined with a heatmap to focus on specific genes instead of their cytobands.
The sample graphs now reside in:

-   Analysis/Statistics/Normals/byBinHeatmap - for Normals

-   Analysis/Statistics/LabLMS/byBinHeatmap - for LabLMS

-   Analysis/Statistics/LM/byBinHeatmap - for LM

It has not been processed for TCGA data.
The graphs have not included methylmaster data due to current constraints.
This section will be updated when it does.
The code for this query and labeling can be found in Visualization.R in the AnalysisCode directory.
