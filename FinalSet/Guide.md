---
editor_options: 
  markdown: 
    wrap: 72
---

# Final Dataset/Analysis Guide

Hi there, this is an informal guide to the results/analysis guide, to
enable effective navigation of the FinalSet folder. It is pretty simple.
All folders are the same names as the titles, but may be missing some
spaces.

### Final Analysis Code

This is a file that contains all the possible lines you could need to
run the full analysis entirely independently. It moves from one analysis
to the other, and is fully signposted as to what each bit of the file
does. You can continue clicking run, and as long as you have the primary
datasets, and can link them properly to the filestructure reader, you
will be able to run it. The prerequisites to run this type of analysis
is to download the repository.

### Accuracies

This uses the accuracy analytical routine to output accuracies across
LMS and LM, with the raw data being stored in the .csv. This is
processed with both mean and median, so the "median" in the name
indicates that it is processed using the median as the summarization
function, while the ones without that marker have used mean. Metrics.png
shows all the data catagories by bins and tech.

### Gene Concordance

It contains a single image, which attempts to quantify the relationships
in genes commonly observed in LMS, by taking the Pearson's Correlation
Coefficient for genes across the 13 cases given. Each of the genes has
their corresponding .csv file for raw data observation.

### Genome Changed

It is pretty simple, it sorts by modified segments, takes its width,
divides it by the basepair width of the entire hg19 genome, and voilà,
we're done. It is then plotted for each technology across each bin for
each case in the GChangedLMS.png, with the raw data being the
corresponding .csv file.

### Genomic Index

This is given by the formula: $$
GI = CNV^2/Cr
$$ Where Cr is the number of chromosomes that the CNVs are observed on.
This is computed only for LMS, with the CorrectedLMSPlot.png plotting the
Genomic index across each technology, bin and case. The recomputed graph
has been added. The high GI actually seems to be a feature of some of
the cases. 
In the case of LM, sometimes the filter seems to be catching that there are just 
no CNVs larger than a certain size, which is expected in a more focal CNV
environment. The associated files are LMGenomicIndex.csv, and the plot is
LMGenomicIndexPlot.png.

### Individual Case Graphs

So the folders are sorted by the sample type, either normal, LMS or LM,
and then below that, it is the bin size. The bin size studied is based
on the given sizes being studied, with the two larger ones being
represented in scientific notation. Under those will be pngs maked with
their STT code. Eahc of those contain a genomic displacement plot (which
is to be fixed to make it publish-worthy), and a gene concordance for
that specific case.
