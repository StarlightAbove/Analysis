library(ggplot2)
library(dplyr)

plot_cnv_segments <- function(df) {
  
  df <- df %>%
    mutate(
      type_group = ifelse(type == "SNP", "SNP", "non-SNP")
    ) %>%
    arrange(chrom, loc.start)
  
  # Calculate chromosome cumulative positions
  chr_lengths <- df %>%
    group_by(chrom) %>%
    summarize(chr_len = max(loc.end), .groups = "drop") %>%
    arrange(chrom) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
    mutate(chr_mid = chr_start + chr_len / 2)
  # Join to get cumulative start and end positions
  df <- df %>%
    left_join(chr_lengths, by = "chrom") %>%
    mutate(
      start_cum = loc.start + chr_start,
      end_cum = loc.end + chr_start
    ) %>%mutate(type = ifelse(type == "SNP", "SNP Array", type)) %>% mutate(type = as.factor(type))
  
  # Vertical chromosome boundaries
  chr_boundaries <- chr_lengths %>%
    mutate(x = chr_start) %>%
    select(chrom, x)
  
  x_breaks <- chr_lengths$chr_mid
  x_labels <- paste0("chr", chr_lengths$chrom)
  
  print(df)
  
  # Plot
  p <- ggplot(df, aes(x = start_cum, xend = end_cum, y = seg.mean, yend = seg.mean)) +
    geom_segment(aes(color = type), size = 0.7, alpha = 0.8) +
    geom_vline(data = chr_boundaries, aes(xintercept = x), color = "grey70", linetype = "dashed") +
    # scale_color_manual(values = c("Amplification" = "red", "Deletion" = "blue", "Normal" = "black")) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = c("SeSAMe" = "red", "MethylMaster" = "blue", "Conumee" = "green", "SNP Array" = "black")) +
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
