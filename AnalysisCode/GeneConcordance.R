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
ap <- rbind(ap10000, ap50000, ap1e05, ap1e06) %>% dplyr::select(-c("chrom", "loc.start", "loc.end"))
ap <- split(ap, ap$Gene)

for (i in seq_along(ap)) {
  file_name <- paste0(names(ap)[i], ".csv")
  write.csv(ap[[i]], file = paste0("~/Work/Analysis/Statistics/LabLMS/GeneConcordance/rawData/", file_name), row.names = FALSE)
}
