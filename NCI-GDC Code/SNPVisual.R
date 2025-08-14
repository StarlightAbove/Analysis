# Visualizing SNP data. (Mainly for data exploration)
library(ggplot2)
library(dplyr)
library(readr)

sexDet <- function(dataFile) {
  if(dataFile[dataFile$Chromosome == "Y", ]$Segment_Mean < -0.3){
    return("Female")
  } else {
    return("Male")
  }
}

# CNV data has columns: Chromosome, Start, End, Probe_Num, Segment_Mean
amp <- 0.3
del <- -0.3
cnv1 <- read_csv(file = "./cases/TCGA-3B-A9HQ/cnvs/cnvsce80d87f-c14a-46e1-a144-5be94e8b16ba.csv", trim_ws = T, progress = T) 
cnv1 <- cnv1[-c(1)]
cnv1[cnv1$Chromosome == "Y", ]
cnv1 <- cnv1 %>% mutate(cnv_type = case_when(Segment_Mean > amp ~ "Amplification",
                                             Segment_Mean < del ~ "Deletion",
                                             TRUE ~ "Normal"
                                             ),
                        copyNumber = 2 * 2^Segment_Mean, 
                        )
ggplot(cnv1) +
  geom_segment(aes(x = Start/1e6, xend = End/1e6, 
                   y = Segment_Mean, yend = Segment_Mean,
                   color = cnv_type, size = log10(Num_Probes + 1))) +
  geom_hline(yintercept = c(del, amp), 
             linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, color = "black") +
  ggtitle(sexDet(cnv1)) +
  scale_color_manual(values = c("Amplification" = "red", 
                                "Deletion" = "blue", 
                                "Normal" = "gray60")) +
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 6) +
  labs(x = "Position (Mb)", y = "Segment Mean (log2 ratio)",
       color = "CNV Type", size = "Segment Weight") +
  theme_minimal()

# CNV analysis.

