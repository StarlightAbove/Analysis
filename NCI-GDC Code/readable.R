# Processing .txt CNV data to readable .csv for comparison, with internal proper 
# labelling.
# Also processing .txt methylation data.
# Side note: Who decided that .txt files were appropriate for this at the NIH?

library("tidyverse")
basePath <- "./cases/"
cases <- list.files(basePath)
paths <- list.files(path = basePath, pattern = ".txt$", recursive = TRUE)
listf <- paste0(basePath, paths)
for(l in listf){
  b <- read.delim(l)
  pathAdd <- str_split(l, "/")
  pathAdd2 <- paste0(trimws(pathAdd[[1]][1:4]), collapse = "/")
  if(grepl("methylation", pathAdd2, fixed = TRUE)){
    write.csv(b, paste0(pathAdd2, "/methylation.csv", collapse = "/"))
  }
}

for(l in listf){
  b <- read.delim(l)
  pathAdd <- str_split(l, "/")
  pathAdd2 <- paste0(trimws(pathAdd[[1]][1:4]), collapse = "/")
  if(grepl("cnvs", pathAdd2, fixed = TRUE)){
    write.csv(b, paste0(pathAdd2, paste0("/cnvs", pathAdd[[1]][5], ".csv", collapse = ""), collapse = "/"))
  }
}

