# Using 4 different software systems, I want to analyse the CNVs from CEL files 
# to establish a highly verified truth dataset.
# DNACopy, QuantiSNP, PennCNV, EnsembleCNV.

# Vignette reproduction.
library(DNAcopy)
library(affy)
data(coriell)
CNA.object <- CNA(cbind(coriell$Coriell.05296),
                  coriell$Chromosome,coriell$Position,
                  data.type="logratio",sampleid="c05296")
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
plot(segment.smoothed.CNA.object, plot.type="s")

# Analysis pipeline.
# Takes log2ratios, with a chromosome, position and a sample id.
data <- affyio::read.celfile(filename = "./LabData/LM_SNP_EPIC_array_data/SNP_array_data_LM/CEL_files/CEL/a520780-00-949463-073116-4238195-70390.CEL")
