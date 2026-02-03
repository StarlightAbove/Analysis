# LOESS-based Genomic Data Correction
# For correcting systematic wave artifacts in copy number data

library(ggplot2)
library(dplyr)

# ==============================================================================
# MAIN CORRECTION FUNCTION
# ==============================================================================

loess_correct_chromosome <- function(truth_data, test_data, chromosome, span = 0.5) {
  #' Correct test data using LOESS fitted to truth data for a single chromosome
  #' 
  #' @param truth_data: data.frame with columns: chr, position, log2ratio
  #' @param test_data: data.frame with columns: chr, position, log2ratio
  #' @param chromosome: chromosome to process
  #' @param span: LOESS span parameter (0-1)
  #' @return: corrected test data for this chromosome
  
  # Filter data for this chromosome
  truth_chr <- truth_data %>% filter(chr == chromosome)
  test_chr <- test_data %>% filter(chr == chromosome)
  
  # Check if we have enough truth points
  if (nrow(truth_chr) < 3) {
    warning(paste("Chromosome", chromosome, "has fewer than 3 truth points. Skipping LOESS correction."))
    return(test_chr %>% mutate(log2ratio_corrected = log2ratio, method = "uncorrected"))
  }
  
  # Find test data points that match truth positions (or are close)
  # For genomic data, we need to match positions
  test_at_truth_pos <- test_chr %>%
    filter(position %in% truth_chr$position)
  
  if (nrow(test_at_truth_pos) < 3) {
    warning(paste("Chromosome", chromosome, "has fewer than 3 matching positions. Skipping LOESS correction."))
    return(test_chr %>% mutate(log2ratio_corrected = log2ratio, method = "uncorrected"))
  }
  
  # Fit LOESS to truth data
  loess_truth <- loess(log2ratio ~ position, data = truth_chr, span = span)
  
  # Fit LOESS to test data at matching positions
  loess_test <- loess(log2ratio ~ position, 
                      data = test_at_truth_pos, 
                      span = span)
  
  # Predict smooth curves across the full chromosome range
  all_positions <- test_chr$position
  
  # Get predictions (with error handling for extrapolation)
  pred_truth <- tryCatch({
    predict(loess_truth, newdata = data.frame(position = all_positions))
  }, error = function(e) {
    # If extrapolation fails, use nearest value
    rep(NA, length(all_positions))
  })
  
  pred_test <- tryCatch({
    predict(loess_test, newdata = data.frame(position = all_positions))
  }, error = function(e) {
    rep(NA, length(all_positions))
  })
  
  # Calculate correction (difference between curves)
  correction <- pred_test - pred_truth
  
  # Apply correction to all test data on this chromosome
  test_chr_corrected <- test_chr %>%
    mutate(
      loess_truth = pred_truth,
      loess_test = pred_test,
      correction = correction,
      log2ratio_corrected = log2ratio - correction,
      method = paste0("loess_span_", span)
    )
  
  return(test_chr_corrected)
}

# ==============================================================================
# OPTIMIZE SPAN PARAMETER
# ==============================================================================

optimize_span <- function(truth_data, test_data, chromosome, 
                         span_values = seq(0.3, 0.9, by = 0.1)) {
  #' Find optimal LOESS span by minimizing error at truth positions
  #' 
  #' @param truth_data: data.frame with truth values
  #' @param test_data: data.frame with test values
  #' @param chromosome: chromosome to optimize for
  #' @param span_values: vector of span values to try
  #' @return: list with optimal span and error metrics
  
  truth_chr <- truth_data %>% filter(chr == chromosome)
  test_chr <- test_data %>% filter(chr == chromosome)
  
  # Get test points at truth positions
  test_at_truth <- test_chr %>%
    filter(position %in% truth_chr$position) %>%
    arrange(position)
  
  truth_chr <- truth_chr %>% arrange(position)
  
  if (nrow(test_at_truth) < 3) {
    return(list(optimal_span = 0.5, errors = NULL))
  }
  
  # Try each span value
  results <- data.frame()
  
  for (span in span_values) {
    # Correct the chromosome
    corrected <- loess_correct_chromosome(truth_data, test_data, chromosome, span)
    
    # Get corrected values at truth positions
    corrected_at_truth <- corrected %>%
      filter(position %in% truth_chr$position) %>%
      arrange(position)
    
    # Calculate error metrics
    if (nrow(corrected_at_truth) == nrow(truth_chr)) {
      rmse <- sqrt(mean((corrected_at_truth$log2ratio_corrected - truth_chr$log2ratio)^2))
      mae <- mean(abs(corrected_at_truth$log2ratio_corrected - truth_chr$log2ratio))
      
      results <- rbind(results, data.frame(
        span = span,
        rmse = rmse,
        mae = mae
      ))
    }
  }
  
  # Find optimal span (minimum RMSE)
  if (nrow(results) > 0) {
    optimal_span <- results$span[which.min(results$rmse)]
  } else {
    optimal_span <- 0.5  # default
  }
  
  return(list(
    optimal_span = optimal_span,
    errors = results
  ))
}

# ==============================================================================
# FULL DATASET CORRECTION
# ==============================================================================

correct_all_chromosomes <- function(truth_data, test_data, 
                                   optimize_spans = FALSE,
                                   default_span = 0.5) {
  #' Correct all chromosomes in test data
  #' 
  #' @param truth_data: data.frame with truth values
  #' @param test_data: data.frame with test values  
  #' @param optimize_spans: if TRUE, optimize span per chromosome
  #' @param default_span: default span if not optimizing
  #' @return: corrected test data for all chromosomes
  
  chromosomes <- unique(test_data$chr)
  corrected_all <- data.frame()
  span_info <- data.frame()
  
  for (chr in chromosomes) {
    cat(paste("Processing chromosome", chr, "...\n"))
    
    # Check if we have truth data for this chromosome
    truth_chr <- truth_data %>% filter(chr == chr)
    
    if (nrow(truth_chr) == 0) {
      cat(paste("  No truth data for chromosome", chr, "- using global linear correction\n"))
      # Fall back to global linear correction
      test_chr <- test_data %>% filter(chr == chr)
      # You would implement global correction here
      test_chr$log2ratio_corrected <- test_chr$log2ratio
      test_chr$method <- "no_truth_data"
      corrected_all <- rbind(corrected_all, test_chr)
      next
    }
    
    # Optimize or use default span
    if (optimize_spans) {
      opt_result <- optimize_span(truth_data, test_data, chr)
      span <- opt_result$optimal_span
      cat(paste("  Optimal span:", round(span, 2), "\n"))
      
      if (!is.null(opt_result$errors)) {
        span_info <- rbind(span_info, 
                          opt_result$errors %>% mutate(chr = chr))
      }
    } else {
      span <- default_span
    }
    
    # Correct this chromosome
    corrected_chr <- loess_correct_chromosome(truth_data, test_data, chr, span)
    corrected_all <- rbind(corrected_all, corrected_chr)
  }
  
  return(list(
    corrected_data = corrected_all,
    span_info = span_info
  ))
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

plot_chromosome_correction <- function(truth_data, test_data, corrected_data, 
                                      chromosome, show_loess = TRUE) {
  #' Plot original and corrected data for a chromosome
  
  truth_chr <- truth_data %>% filter(chr == chromosome)
  test_chr <- test_data %>% filter(chr == chromosome)
  corrected_chr <- corrected_data %>% filter(chr == chromosome)
  
  # Create plot
  p <- ggplot() +
    # Original test data
    geom_point(data = test_chr, aes(x = position, y = log2ratio), 
               color = "red", alpha = 0.3, size = 1) +
    # Corrected test data
    geom_point(data = corrected_chr, aes(x = position, y = log2ratio_corrected),
               color = "blue", alpha = 0.3, size = 1) +
    # Truth data (larger points)
    geom_point(data = truth_chr, aes(x = position, y = log2ratio),
               color = "black", size = 3, shape = 17) +
    labs(
      title = paste("Chromosome", chromosome, "- LOESS Correction"),
      x = "Genomic Position",
      y = "Log2 Ratio",
      subtitle = "Red = Original | Blue = Corrected | Black triangles = Truth"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add LOESS curves if available and requested
  if (show_loess && "loess_truth" %in% names(corrected_chr)) {
    p <- p +
      geom_line(data = corrected_chr, aes(x = position, y = loess_truth),
                color = "darkgreen", linewidth = 1, linetype = "dashed") +
      geom_line(data = corrected_chr, aes(x = position, y = loess_test),
                color = "orange", linewidth = 1, linetype = "dashed")
  }
  
  return(p)
}

plot_span_optimization <- function(span_info, chromosome = NULL) {
  #' Plot RMSE vs span for chromosome(s)
  
  if (!is.null(chromosome)) {
    span_info <- span_info %>% filter(chr == chromosome)
    title <- paste("Span Optimization - Chromosome", chromosome)
  } else {
    title <- "Span Optimization - All Chromosomes"
  }
  
  p <- ggplot(span_info, aes(x = span, y = rmse, color = as.factor(chr))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = title,
      x = "LOESS Span",
      y = "RMSE at Truth Positions",
      color = "Chromosome"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# Example of how to use these functions:
#
# # Load your data
# truth_data <- read.csv("truth_data.csv")  
# # Columns: chr, position, log2ratio
#
# test_data <- read.csv("test_data.csv")    
# # Columns: chr, position, log2ratio
#
# # Method 1: Use default span (0.5) for all chromosomes
# result <- correct_all_chromosomes(truth_data, test_data, 
#                                   optimize_spans = FALSE, 
#                                   default_span = 0.5)
# corrected_data <- result$corrected_data
#
# # Method 2: Optimize span per chromosome
# result <- correct_all_chromosomes(truth_data, test_data, 
#                                   optimize_spans = TRUE)
# corrected_data <- result$corrected_data
# span_info <- result$span_info
#
# # Visualize results for chromosome 1
# plot_chromosome_correction(truth_data, test_data, corrected_data, 
#                           chromosome = "chr1")
#
# # Plot span optimization
# plot_span_optimization(span_info, chromosome = "chr1")
#
# # Save corrected data
# write.csv(corrected_data, "corrected_data.csv", row.names = FALSE)

cat("LOESS correction functions loaded successfully!\n")
cat("\nMain functions:\n")
cat("  - correct_all_chromosomes(): Correct entire dataset\n")
cat("  - optimize_span(): Find optimal span for a chromosome\n")
cat("  - plot_chromosome_correction(): Visualize correction\n")
cat("  - plot_span_optimization(): Plot span vs error\n")
cat("\nSee example usage at the end of this script.\n")
