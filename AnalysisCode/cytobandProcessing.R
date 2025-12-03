reference <- data.frame(
  Cytoband = c("8q24.21", "17p11.2", "9q12", "9q34.3", "10q23.31", "13q14.2", 
               "17p13.1"),
  Gene = c("MYC", "MYOCD", "CCNE1", "CDKN2A", "PTEN", "RB1", "TP53"),
  CNVStatus = c("Deletion", "Deletion","Deletion","Amplification", 
                     "Amplification", "Amplification", "Amplification"),
  seg.mean = c(-0.2, -0.2, -0.2, 0.2, 0.2, 0.2, 0.2),
  type = c("Gene")
) 
reference <- reference %>% rowwise() %>% 
  mutate(loc.start = as.numeric(convert_band_to_genomic(Cytoband)[["start"]][["loc.start"]]),
          loc.end = as.numeric(convert_band_to_genomic(Cytoband)[["end"]][["loc.end"]]),
          chrom = as.numeric(str_remove_all(convert_band_to_genomic(Cytoband)[["chromosome"]][["Chromosome"]], 
                                            "chr")),
         type = as.factor(type)) %>% 
  ungroup() %>% select(c("CNVStatus", "loc.start", "loc.end", "seg.mean", "chrom", "type"))

data(cytobandLocations) # This loads the 'cytobandLocations' data frame

# View a sample of the data structure (optional)
# head(cytobandLocations)

# 2. Define a function to convert a given band name
convert_band_to_genomic <- function(band_name) {
  # Ensure the input name is in the correct format (e.g., "Xq28")
  # The 'cytobandLocations' uses names like "Xq28" as rownames
  
  if (!(band_name %in% rownames(cytobandLocations))) {
    return(paste("Band", band_name, "not found in the reference data."))
  }
  
  # Retrieve the row corresponding to the band name
  band_info <- cytobandLocations[band_name, ]
  
  # Extract coordinates
  chromosome <- band_info["Chromosome"]
  start_pos <- band_info["loc.start"]
  end_pos <- band_info["loc.end"]
  
  return(list(
    chromosome = chromosome,
    start = start_pos,
    end = end_pos,
    genome_build = "GRCh38/hg19" # Data is compatible with both
  ))
}
s <- convert_band_to_genomic("8q24.21")
