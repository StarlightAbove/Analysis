library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(cowplot)
library(gridGraphics)
library(gridExtra)
library(ggplotify)
library(AnnotationHub)
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("rentrez")
library("pheatmap")
library("GenomicRanges")
library(RColorBrewer)
hs <- org.Hs.eg.db

geneAnno <- function(Gene, db = NULL){
  reference <- geneGen(Gene = Gene, db = db)
  print(reference[[1]])
  pheatmap_ggplot <- as.ggplot(reference[[1]]$gtable)
  reference2 <- reference[[2]]
  db <- unique(rbind(db, reference2))
  img2 <- plot_cnv_segments(db)
  return(list((img2 + pheatmap_ggplot), db))
}
plot_cnv_segments <- function(df, anno = NULL) {
  
  df <- df %>%
    mutate(
      type_group = ifelse(type == "SNP", "SNP", "non-SNP")
    ) %>%
    arrange(chrom, loc.start)
  
  # Calculate chromosome cumulative positions
  chr_lengths <- df %>%
    dplyr::group_by(chrom) %>%
    dplyr::summarize(chr_len = max(loc.end), .groups = "drop") %>%
    arrange(chrom) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
    mutate(chr_mid = chr_start + chr_len / 2)
  print(chr_lengths)
  # Join to get cumulative start and end positions
  df <- df %>%
    left_join(chr_lengths, by = "chrom") %>%
    mutate(
      start_cum = loc.start + chr_start,
      end_cum = loc.end + chr_start
    ) %>%mutate(type = ifelse(type == "SNP", "SNP Array", type)) %>% 
    mutate(type = as.factor(type))
  
  # Vertical chromosome boundaries
  chr_boundaries <- chr_lengths %>%
    mutate(x = chr_start) %>%
    dplyr::select(chrom, x)
  
  x_breaks <- chr_lengths$chr_mid
  x_labels <- paste0("chr", chr_lengths$chrom)
  x_labels <- c(x_labels)
  print(x_labels)
  print(df)
  df$Gene[df$Gene == "0"] <- NA
  
  # Plot
  p <- ggplot(df, aes(x = start_cum, xend = end_cum, y = seg.mean, yend = seg.mean, label = Gene)) +
    geom_segment(aes(color = type), size = 0.7, alpha = 0.8) +
    geom_point(data = subset(df, type == "Gene_Conumee"), aes(x = start_cum, y = seg.mean, label = Gene), color = "green", size = 3) + 
    geom_point(data = subset(df, type == "Gene_SNP"), aes(x = start_cum, y = seg.mean, label = Gene), color = "gray40", size = 3) + 
    geom_point(data = subset(df, type == "Gene_SeSAMe"), aes(x = start_cum, y = seg.mean, label = Gene), color = "red", size = 3) + 
    geom_text_repel(
      box.padding = unit(0.35, "lines"), # Adjust padding around the label
      point.padding = unit(0.3, "lines"), # Adjust padding around the point
      segment.color = 'grey' # Color of the connecting lines
    ) +
    #geom_text(data = subset(df, type == "Gene"), mapping = aes(x = start_cum,
    # y = 0,
    # label = Gene,
    # hjust = -1,
    # vjust = -1)) + 
    geom_vline(data = chr_boundaries, aes(xintercept = x), color = "grey70", linetype = "dashed") +
    # scale_color_manual(values = c("Amplification" = "red", "Deletion" = "blue", "Normal" = "black")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("SeSAMe" = "red", "MethylMaster" = "blue", 
                                  "Conumee" = "green", "SNP Array" = "black", "Gene" = "orange")) +
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
  
  return(p)
}
# Dependency for geneAnno
Gene <- c("MYC", "MYOCD", "CCNE1", "CDKN2A", "PTEN", "RB1", "TP53") 
geneGen <- function(Gene, db = NULL){
  entrez_idsDB <- mapIds(org.Hs.eg.db, 
                         keys=Gene, 
                         column="ENTREZID", 
                         keytype="SYMBOL",
                         multiVals="first")
  chr_locations <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = entrez_idsDB,
                                         columns = c("CHR", "CHRLOC", "CHRLOCEND"),
                                         keytype = "ENTREZID")
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  target_genes_ranges <- gene_ranges[gene_ranges$gene_id %in% entrez_ids]
  coords_df <- as.data.frame(target_genes_ranges)
  coords_df$entrez_id <- target_genes_ranges$gene_id
  coords_final <- coords_df[, c("entrez_id", "seqnames", "start", "end", 
                                "strand", "width")] %>% 
    dplyr::mutate(chrom=seqnames) %>% 
    dplyr::mutate(chrom = as.numeric(gsub("chr","",chrom)), 
                  Gene = rownames(coords_final)) %>% 
    dplyr::select(-c(seqnames))
  rownames(coords_final) <- Gene
  if(!is.null(db)){
    print("is here")
    dt <- db %>% dplyr::filter((Gene == "0")) %>% 
      dplyr::mutate(Gene = NA) 
    gcoords <- coords_final %>% 
      dplyr::mutate(loc.start = as.numeric(start), loc.end = as.numeric(end)) %>% 
      dplyr::mutate(seg.mean = 0, CNVStatus = "Normal", type = "Gene") %>%
      dplyr::mutate(Gene = Gene) %>% 
      dplyr::select(-c("strand", "width", "entrez_id", "start", "end"))
    grC <- makeGRangesFromDataFrame(gcoords, 
                                    start.field = "loc.start",
                                    end.field = "loc.end",
                                    keep.extra.columns = TRUE)
    grDt <- makeGRangesFromDataFrame(dt,
                                     start.field = "loc.start",
                                     end.field = "loc.end",
                                     keep.extra.columns = TRUE)
    overlaps <- findOverlaps(grDt, grC)
    overlapping_regions <- pintersect(grDt[queryHits(overlaps)], grC[subjectHits(overlaps)])
    overlapping_regions <- as.data.frame(overlapping_regions) %>% 
      dplyr::select(-c("strand", "width", "hit", "Gene")) %>%
      dplyr::rename(loc.start = start, loc.end = end) %>%
      dplyr::mutate(chrom = as.numeric(seqnames)) %>%
      dplyr::select(-c("seqnames"))
    overlapping_regions <- dplyr::left_join(overlapping_regions, coords_final, 
                                            by=c("chrom"="chrom",
                                                 "loc.start" = "start", 
                                                 "loc.end" = "end")) %>%
      dplyr::select(-c("entrez_id", "strand", "width")) 
    overlapping_regions <- unique(overlapping_regions[c(6, 1, 2, 3, 4, 5, 7)])
    overlapping_regions$type <- ifelse(
      overlapping_regions$type == "Conumee", "Gene_Conumee",
      ifelse(overlapping_regions$type == "SeSAMe", "Gene_SeSAMe", 
             ifelse(overlapping_regions$type == "MethylMasteR", "Gene_MMasteR",
                    ifelse(overlapping_regions$type == "SNP", "Gene_SNP", "Unknown")))
    )
    print(overlapping_regions)
    # Developing visualization
    mat <- matrix(NA, nrow = 3, ncol = 7)
    colnames(mat) <- unique(overlapping_regions$Gene)
    rownames(mat) <- unique(overlapping_regions$type)
    for (i in 1:nrow(overlapping_regions)) {
      row_index <- match(overlapping_regions$type[i], rownames(mat))
      col_index <- match(overlapping_regions$Gene[i], colnames(mat))
      mat[row_index, col_index] <- overlapping_regions$seg.mean[i]
    }
    mat[is.na(mat)] <- 0
    print(mat)
    colors <- rev(brewer.pal(n = 7, name = "RdBu"))
    labels_matrix <- matrix(sprintf("%.3f", mat), 
                            nrow = nrow(mat), 
                            ncol = ncol(mat))
    ph <- pheatmap(
      mat,
      color = colors,
      scale = "row", # Scales the values in each row/column/none to a z-score
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = TRUE, # Set to TRUE if you have few genes
      display_numbers = labels_matrix,
      main = "Gene Expression Heatmap (Log2 Ratio / Z-score)"
    )
    list(ph, overlapping_regions)
    
  } else {
    coords_final
  }
}


# Sample

geneGen(Gene = Gene)

# ‘cnAnalysis450k’ are not available for package ‘MethylMasteR’

