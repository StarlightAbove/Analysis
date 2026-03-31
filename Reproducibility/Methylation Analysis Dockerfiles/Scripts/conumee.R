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