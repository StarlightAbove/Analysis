# Necessary function from EaCoN for diploid recentering.
chromobjector <- function(BSg = NULL) {
  if (is.null(BSg)) stop("NULL object !", call. = FALSE)
  # chromobj <- list(species = GenomeInfoDb::organism(BSg), genomebuild = BSgenome::providerVersion(BSg))
  chromobj <- list(species = GenomeInfoDb::organism(BSg), genomebuild = metadata(BSg)$genome)
  chromdf <- data.frame(chrom = BSgenome::seqnames(BSg), chrN = seq_along(BSgenome::seqnames(BSg)), chr.length = GenomeInfoDb::seqlengths(BSg), stringsAsFactors = FALSE)
  chromdf$chr.length.sum <- cumsum(as.numeric(chromdf$chr.length))
  chromdf$chr.length.toadd <- c(0, chromdf$chr.length.sum[-nrow(chromdf)])
  chromdf$mid.chr <- round(diff(c(0, chromdf$chr.length.sum)) /2)
  chromdf$mid.chr.geno <- chromdf$mid.chr + chromdf$chr.length.toadd
  chromobj$chromosomes <- chromdf
  rm(chromdf)
  chromobj$chrom2chr <- sapply(chromobj$chromosomes$chrom, function(k) { chromobj$chromosomes$chrN[chromobj$chromosomes$chrom == k]}, simplify = FALSE)
  chromobj$chr2chrom <- sapply(chromobj$chromosomes$chrN, function(k) { chromobj$chromosomes$chrom[chromobj$chromosomes$chrN == k]}, simplify = FALSE)
  names(chromobj$chr2chrom) <- chromobj$chromosomes$chrN
  chromobj$genome.length <- sum(as.numeric(chromobj$chromosomes$chr.length), na.rm = TRUE)
  return(chromobj)
}
