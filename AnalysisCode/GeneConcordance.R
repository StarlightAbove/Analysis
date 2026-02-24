# Generate a table for each gene
# Column = tech + SNP
# Row = cases
# Dataset = log2ratio
library(ggthemes)
# LMS
Gene <- c("MYC", "MYOCD", "CCNE1", "CDKN2A", "PTEN", "RB1", "TP53") 
stts <- c(9202, 9203, 9327, 9328, 9337, 9338, 9350, 9353, 9354, 9355, 9356, 9357, 9358)
# Problem candidates: 9327 @ 1e+05 and 1e+06
bins <- c(10000, 50000, 1e+05, 1e+06)

NoGraphGeneGen <- function(Gene, db = NULL, case){
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
             ifelse(overlapping_regions$type == "MethylMaster", "Gene_MMasteR",
                    ifelse(overlapping_regions$type == "SNP", "Gene_SNP", "Unknown")))
    )
    overlapping_regions <- overlapping_regions %>% dplyr::mutate(case = case)
    overlapping_regions
    
  } else {
    coords_final
  }
}
NoGraphGeneGen <- function(Gene, db = NULL, case){
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
                  Gene = Gene) %>% 
    dplyr::select(-c(seqnames))
  rownames(coords_final) <- Gene
  if(!is.null(db)){
    print("is here")
    dt <- db %>% dplyr::filter((Gene == "0")) %>% 
      dplyr::mutate(Gene = NA)
    print(dt %>% dplyr::filter(type == "MethylMaster"))
    gcoords <- coords_final %>% 
      dplyr::mutate(loc.start = as.numeric(start), loc.end = as.numeric(end)) %>% 
      dplyr::mutate(seg.mean = 0, CNVStatus = "Normal", type = "Gene") %>%
      dplyr::mutate(Gene = Gene) %>% 
      dplyr::select(-c("strand", "width", "entrez_id", "start", "end"))
    print(gcoords)
    grC <- makeGRangesFromDataFrame(gcoords, 
                                    start.field = "loc.start",
                                    end.field = "loc.end",
                                    keep.extra.columns = TRUE)
    grDt <- makeGRangesFromDataFrame(dt,
                                     start.field = "loc.start",
                                     end.field = "loc.end",
                                     keep.extra.columns = TRUE)
    
    print("grDt chrom 19 MethylMasteR rows:")
    print(as.data.frame(grDt) %>% dplyr::filter(seqnames == 19, type == "MethylMaster"))
    
    overlaps <- findOverlaps(grDt, grC)
    
    print("Overlap hits for MethylMasteR chrom 19:")
    print(as.data.frame(grDt[queryHits(overlaps)]) %>% dplyr::filter(seqnames == 19, type == "MethylMaster"))
    
    
    # Extract intersection and bind gene annotation immediately before any row changes
    intersection_df <- as.data.frame(pintersect(
      grDt[queryHits(overlaps)], 
      grC[subjectHits(overlaps)]
    )) %>%
      dplyr::select(-c("strand", "width", "hit", "Gene")) %>%
      dplyr::rename(loc.start = start, loc.end = end) %>%
      dplyr::mutate(chrom = as.numeric(seqnames)) %>%
      dplyr::select(-c("seqnames"))
    
    gene_annotations <- coords_final[subjectHits(overlaps), "Gene", drop = FALSE]
    overlapping_regions <- cbind(intersection_df, Gene = gene_annotations$Gene)
    
    # Find genes with no overlap for each type
    all_types <- unique(dt$type)
    missing_rows <- list()
    
    for(t in all_types){
      dt_type <- dt %>% dplyr::filter(type == t)
      grDt_type <- makeGRangesFromDataFrame(dt_type,
                                            start.field = "loc.start",
                                            end.field = "loc.end",
                                            keep.extra.columns = TRUE)
      matched_genes <- unique(overlapping_regions$Gene[overlapping_regions$type == t])
      missing_genes <- setdiff(rownames(coords_final), matched_genes)
      
      for(g in missing_genes){
        gene_row <- coords_final[g, ]
        gene_mid <- (gene_row$start + gene_row$end) / 2
        
        # Find flanking segments on same chromosome
        dt_chr <- dt_type %>% dplyr::filter(chrom == gene_row$chrom)
        
        if(nrow(dt_chr) == 0) next
        
        # Left flank: segments ending before gene midpoint
        left <- dt_chr %>% dplyr::filter(loc.end < gene_mid) %>% 
          dplyr::slice_max(loc.end, n = 1)
        # Right flank: segments starting after gene midpoint  
        right <- dt_chr %>% dplyr::filter(loc.start > gene_mid) %>% 
          dplyr::slice_min(loc.start, n = 1)
        
        if(nrow(left) == 0 || nrow(right) == 0) next
        
        avg_seg <- mean(c(left$seg.mean, right$seg.mean))
        avg_cnv <- ifelse(left$CNVStatus == right$CNVStatus, 
                          left$CNVStatus, "Normal")
        
        missing_rows[[length(missing_rows) + 1]] <- data.frame(
          chrom = gene_row$chrom,
          loc.start = gene_row$start,
          loc.end = gene_row$end,
          seg.mean = avg_seg,
          type = t,
          CNVStatus = avg_cnv,
          Gene = g,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if(length(missing_rows) > 0){
      overlapping_regions <- rbind(overlapping_regions, 
                                   dplyr::bind_rows(missing_rows))
    }
    
    print("intersection_df chrom 19 MethylMasteR rows:")
    print(intersection_df %>% dplyr::filter(chrom == 19, type == "MethylMaster"))
    
    print("overlapping_regions before unique, chrom 19 MethylMasteR:")
    print(overlapping_regions %>% dplyr::filter(chrom == 19, type == "MethylMaster"))
    
    print(overlapping_regions %>% dplyr::filter(Gene == "RB1"))
    overlapping_regions <- unique(overlapping_regions[c(6, 1, 2, 3, 4, 5, 7)])
    overlapping_regions$type <- ifelse(
      overlapping_regions$type == "Conumee", "Gene_Conumee",
      ifelse(overlapping_regions$type == "SeSAMe", "Gene_SeSAMe", 
             ifelse(overlapping_regions$type == "MethylMaster", "Gene_MMasteR",
                    ifelse(overlapping_regions$type == "SNP", "Gene_SNP", "Unknown")))
    )
    overlapping_regions <- overlapping_regions %>% dplyr::mutate(case = case)
    overlapping_regions
    
  } else {
    coords_final
  }
}
NoGraphGeneGen(Gene = Gene, case = 9357)

ap10000 <- NULL
ap50000 <- NULL
ap1e05 <- NULL
ap1e06 <- NULL

NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(9357, "MethylMaster", 10000), 
                                       labLMSProc(9357, "Conumee", 10000), 
                                       labLMSProc(9357, "Sesame", 10000)), 
               case = 9357) %>% dplyr::filter(Gene == "RB1")


for(stt in stts){
  ap10000 <- rbind(ap10000, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 10000), 
                                                    labLMSProc(stt, "Conumee", 10000), 
                                                    labLMSProc(stt, "Sesame", 10000)), 
                            case = stt))
  ap50000 <- rbind(ap50000, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 50000), 
                                                    labLMSProc(stt, "Conumee", 50000), 
                                                    labLMSProc(stt, "Sesame", 50000)), 
                            case = stt))
  ap1e05 <- rbind(ap1e05, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 1e+05), 
                                                   labLMSProc(stt, "Conumee", 1e+05), 
                                                   labLMSProc(stt, "Sesame", 1e+05)), 
                           case = stt))
  ap1e06 <- rbind(ap1e06, NoGraphGeneGen(Gene = Gene, db = rbind(labLMSProc(stt, "MethylMaster", 1e+06), 
                                                   labLMSProc(stt, "Conumee", 1e+06), 
                                                   labLMSProc(stt, "Sesame", 1e+06)), 
                           case = stt))
}
ap10000 %>% filter(Gene == "RB1", case == 9356)
techs <- c("Gene_MMasteR", "Gene_SNP", "Gene_Conumee", "Gene_SeSAMe")

ap10000 <- ap10000 %>% dplyr::rename(log2ratio = seg.mean)
ap50000 <- ap50000 %>% dplyr::rename(log2ratio = seg.mean)
ap1e05 <- ap1e05 %>% dplyr::rename(log2ratio = seg.mean)
ap1e06 <- ap1e06 %>% dplyr::rename(log2ratio = seg.mean)

ggplot(ap10000, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 10000") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(ap50000, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 50000") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(ap1e05, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 1e+05") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

ggplot(ap1e06, aes(as.factor(case), Gene, fill=log2ratio)) +
  facet_wrap(~type) +
  xlab("Case") +
  ylab("Oncogene") +
  labs(title = "Relationship between log-2 ratio, technology & gene",
       subtitle = "Bin Size: 1e+06") +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()


# Saving raw data
ap10000 <- ap10000 %>% dplyr::mutate(bin = 10000)
ap50000 <- ap50000 %>% dplyr::mutate(bin = 50000)
ap1e05 <- ap1e05 %>% dplyr::mutate(bin = 1e+05)
ap1e06 <- ap1e06 %>% dplyr::mutate(bin = 1e+06)
ap <- rbind(ap10000, ap50000, ap1e05, ap1e06) %>% dplyr::select(-c("chrom")) %>%
  dplyr::mutate(width = loc.end - loc.start)
ap <- split(ap, ap$Gene)

for (i in seq_along(ap)) {
  file_name <- paste0(names(ap)[i], ".csv")
  write.csv(ap[[i]], file = paste0("~/Work/Analysis/Statistics/LabLMS/GeneConcordance/rawData/", file_name), row.names = FALSE)
}
ap[[Gene[1]]] <- split(ap[[Gene[1]]], ap[[Gene[1]]]$bin)
ap[[Gene[2]]] <- split(ap[[Gene[2]]], ap[[Gene[2]]]$bin)
ap[[Gene[3]]] <- split(ap[[Gene[3]]], ap[[Gene[3]]]$bin)
ap[[Gene[4]]] <- split(ap[[Gene[4]]], ap[[Gene[4]]]$bin)
ap[[Gene[5]]] <- split(ap[[Gene[5]]], ap[[Gene[5]]]$bin)
ap[[Gene[6]]] <- split(ap[[Gene[6]]], ap[[Gene[6]]]$bin)
ap[[Gene[7]]] <- split(ap[[Gene[7]]], ap[[Gene[7]]]$bin)

# Looking for concordance
genes <- names(ap)
bins <- c("50000", "10000", "1e+05", "1e+06")
states <- c("Conumee – default bin size",
            "Conumee – 10kb",
            "Conumee – 100kb",
            "Conumee – 1Mb",            
            "Sesame – default bin size",
            "Sesame – 10kb",
            "Sesame – 100kb",
            "Sesame – 1Mb",
            "MethylMasteR – default bin size",
            "MethylMasteR – 10kb",
            "MethylMasteR – 100kb",
            "MethylMasteR – 1Mb"
)

Methyl <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_MMasteR")
Ses <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_SeSAMe")
Con <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_Conumee")
Snp <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_SNP")
Snp <- Snp %>% group_by(case) %>%
  summarize(
    log2ratio = mean(log2ratio)
  )

CCNE1 <- c(1:length(states))
CDKN2A <- c(1:length(states))
MYC <- c(1:length(states))
MYOCD <- c(1:length(states))
PTEN <- c(1:length(states))
RB1 <- c(1:length(states))
TP53 <- c(1:length(states))

output.df <- data.frame(states, CCNE1,
                        CDKN2A,
                        MYC,
                        MYOCD, 
                        PTEN ,
                        RB1 ,
                        TP53)
for(i in seq_along(ap)){
  for(j in bins){
    print(genes[i])
    print(j)
    Methyl <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_MMasteR")
    Ses <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_SeSAMe")
    Con <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_Conumee")
    Snp <- ap[[i]][[j]] %>% dplyr::filter(type == "Gene_SNP")
    Snp <- Snp %>% group_by(case) %>%
      summarize(
        log2ratio = mean(log2ratio)
      )
    mCorr <- cor(Methyl$log2ratio, Snp$log2ratio)
    sCorr <- cor(Ses$log2ratio, Snp$log2ratio)
    cCorr <- cor(Con$log2ratio, Snp$log2ratio)
    print(paste0(mCorr, " ", sCorr, " ", cCorr))
    tag <- 0
    if(j == "50000"){
      tag <- 1
    }
    if(j == "10000"){
      tag <- 2
    }
    if(j == "1e+05"){
      tag <- 3
    }
    if(j == "1e+06"){
      tag <- 4
    }

    output.df[tag, i + 1] <- mCorr
    output.df[tag + 4, i + 1] <- sCorr 
    output.df[tag + 8, i + 1] <- cCorr 
  }
}
# Melt to long format
df_long <- melt(output.df, id.vars = "states", variable.name = "gene", value.name = "correlation")

# Keep row order
df_long$state <- factor(df_long$state, levels = rev(data$state))

# Plot
p <- ggplot(df_long, aes(x = gene, y = states, fill = correlation)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#50C878",
    midpoint = 0.5,
    limits = c(min(df_long$correlation), 1),
    name = "Correlation"
  ) +
  labs(
    title = "CNV Method Comparison by Gene",
    x = "Gene",
    y = "Method & Bin Size"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.position = "right"
  )
  

