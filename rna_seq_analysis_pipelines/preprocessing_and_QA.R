# ==============================================================================
# Author: Noa Gilad
# Title: Preprocessing and QA Utilities for Count-Based RNA-seq Analyses
# Description:
#   Helper functions for preprocessing, quality assessment, visualization, and
#   sample aggregation used across DESeq2 and glmmSeq pipelines. Includes low-count
#   filtering, metadata preparation, count and expression QC plots, outlier scoring,
#   sample-to-sample distance heatmaps, per-sample distribution plots, and utilities
#   to merge technical or biological replicates by group.
#
# Notes:
#   - Functions are designed to be project-agnostic.
#   - Several plotting functions assume log2(Expression + 1) matrices as input.
#   - Gene identifiers are assumed to be ENSEMBL where annotation is added.
#
# Dependencies: dplyr, tidyverse, ggplot2, DESeq2, pheatmap, cowplot, org.Hs.eg.db
# ==============================================================================

# This is a file with functions that are necessary for glmm and for deseq2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Load required libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# tidy-verse style data wrangling ------------------------------------------------
library(dplyr)            # %>%, mutate(), arrange(), join helpers, piping
library(readr)            # read_csv(), write_tsv(), fast Delimiter-based I/O
library(tibble)           # tibble(), as_tibble(); keeps data frames tidy
library(tidyverse)

# visualisation helpers ---------------------------------------------------------
library(ggplot2)          # grammar of graphics base layer
library(ggpubr)           # publication-ready wrappers (stat_compare_means, etc.)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' Filter low-count genes from a count matrix
#'
#' Removes genes with low total read support by keeping only rows whose row sum
#' is greater than a user-defined threshold.
#'
#' @param counts A raw counts matrix with genes as rows and samples as columns.
#' @param readsLim Minimum total reads per gene required to keep the gene.
#'
#' @return A filtered counts matrix.
#'
# a function that prepares the count table (filters low counts)
prep_counts <- function(counts, readsLim) {
  return(counts[which(rowSums(counts) > readsLim),])
}

#' Prepare a colData table from a CSV file
#'
#' Reads a sample design file, converts columns to factors, and relevels the
#' specified condition column to a requested reference level.
#'
#' @param file Path to a CSV design file with samples as rownames in the first column.
#' @param condition_col Column name to treat as the primary condition factor.
#' @param ref_level Reference level for the condition factor.
#'
#' @return A metadata table suitable for use as colData.
#'
# a function that prepares the colData
prep_coldata <- function(file, condition_col, ref_level) {
  design_file <- read.csv(file, sep = ",", row.names = 1, stringsAsFactors = FALSE)
  design_file[] <- lapply(design_file, factor)
  design_file[[condition_col]] <- relevel(design_file[[condition_col]], ref = ref_level)
  return(design_file)
}

#' Create a bar plot of total read counts per sample
#'
#' Computes total reads per sample as the column sums of a count table and
#' returns a bar plot useful for basic sequencing depth quality assessment.
#'
#' @param count_table A raw counts matrix with genes as rows and samples as columns.
#' @param title Plot title.
#' @param color Outline color for bars.
#' @param fill Fill color for bars.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#'
#' @return A ggplot object.
#'
# a function that creates a histogram with gene counts per sample for QA
create_count_histogram <- function(count_table, 
                                   title = "Total Read Count per Sample",
                                   color = "steelblue",
                                   fill = "steelblue",
                                   x_label = "Sample",
                                   y_label = "Total Read Count") {
  
  # Calculate column sums (total reads per sample)
  total_counts <- colSums(count_table)
  
  # Convert to data frame for ggplot
  plot_data <- data.frame(
    Sample = factor(names(total_counts), levels = names(total_counts)),
    Count = total_counts
  )
  
  p <- ggplot(plot_data, aes(x = Sample, y = Count)) +
    geom_bar(stat = "identity", color = color, fill = fill) +
    labs(title = title,
         x = x_label,
         y = y_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


#' Plot pairwise sample expression scatter with 1:1 reference line
#'
#' Creates a scatter plot comparing two samples across genes using a log-scale
#' expression matrix, including a 1:1 line and equal axis scaling for visual QC.
#'
#' @param log_counts A matrix or data.frame of log-scale expression values with samples as columns.
#' @param sample_x Column name for the x-axis sample.
#' @param sample_y Column name for the y-axis sample.
#'
#' @return A ggplot object.
#'
plot_pairwise_expression <- function(log_counts, sample_x, sample_y) {
  df <- data.frame(
    x = log_counts[[sample_x]],
    y = log_counts[[sample_y]]
  )
  
  ggplot(df, aes(x = x, y = y)) +
    # semi-transparent dark points
    geom_point(color = "darkgrey", alpha = 0.4, size = 0.6) +
    # solid black 1:1 line
    geom_abline(slope = 1, intercept = 0, color = "black", size = 0.5) +
    # equal scales on both axes
    coord_equal() +
    # no grid lines, minimal axes
    theme_classic(base_size = 14) +
    theme(
      panel.border     = element_rect(fill = NA, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = sample_x,
      y = sample_y
    )
}

#' Generate and save pairwise expression plots for many sample comparisons
#'
#' For each focus sample, creates pairwise scatter plots versus other samples,
#' groups them into grids, and saves each grid as a PNG in an output directory.
#' Can optionally limit to a subset of samples defined by a metadata grouping.
#'
#' @param log_counts A matrix or data.frame of log-scale expression values with samples as columns.
#' @param coldata Sample metadata table with rownames matching sample names.
#' @param all Logical; if TRUE uses all samples, otherwise subsets by group_col and group_value.
#' @param group_col Column name in coldata used for subsetting samples when all is FALSE.
#' @param group_value Value in group_col used to select samples when all is FALSE.
#' @param out_dir Output directory for saved PNG files.
#' @param max_plots_per_grid Maximum number of pairwise plots per saved grid image.
#'
#' @return No return value. Side effect is saving PNG files to out_dir.
#'
plot_all_pairs <- function(log_counts,
                           coldata,
                           all = TRUE,
                           group_col = NULL,
                           group_value = NULL,
                           out_dir,
                           max_plots_per_grid = 3) {
  # create output directory
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # decide which samples to include
  universe <- if (all || is.null(group_col)) {
    colnames(log_counts)
  } else {
    intersect(
      rownames(coldata)[coldata[[group_col]] == group_value],
      colnames(log_counts)
    )
  }
  
  # for each “focus” sample
  for (s in universe) {
    others <- setdiff(universe, s)
    # split others into chunks of size <= max_plots_per_grid
    chunks <- split(others, ceiling(seq_along(others) / max_plots_per_grid))
    
    # for each chunk, build and save a grid
    for (i in seq_along(chunks)) {
      chunk <- chunks[[i]]
      # make plots for this chunk
      plots <- lapply(chunk, function(o) {
        plot_pairwise_expression(log_counts, s, o)
      })
      names(plots) <- vapply(chunk, function(o) paste(s, o, sep = "_vs_"), character(1))
      
      # assemble with cowplot
      grid_body <- plot_grid(plotlist = plots, ncol = length(plots))
      title     <- ggdraw() +
        draw_label(paste("Comparisons:", s, "vs", paste(chunk, collapse = ", ")),
                   fontface = "bold", x = 0, hjust = 0)
      full_plot <- plot_grid(title, grid_body, ncol = 1, rel_heights = c(0.1, 1))
      
      # save to PNG
      fname <- file.path(out_dir, sprintf("%s_chunk%02d.png", s, i))
      ggsave(fname, full_plot,
             width  = 3 * 3,      
             height = 3,           
             dpi    = 150,
             bg     = "white") 
    }
  }
}


#' Compute average Pearson correlation per sample for outlier detection
#'
#' Computes a pairwise Pearson correlation matrix between samples and returns
#' the average correlation per sample, sorted from lowest to highest.
#'
#' @param log_counts A matrix or data.frame of log-scale expression values with samples as columns.
#' @param coldata Sample metadata table with rownames matching sample names.
#' @param group_col Column name used to subset samples when all is FALSE.
#' @param group_value Value in group_col used to select samples when all is FALSE.
#' @param all Logical; if TRUE uses all samples, otherwise subsets by group_col and group_value.
#'
#' @return A named numeric vector of average correlations sorted ascending.
#'
get_correlation_scores <- function(log_counts, coldata, group_col, group_value, all = FALSE) {
  if (all) {
    # Use all samples
    sub_counts <- log_counts
  } else {
    # Subset to samples from the specified group
    group_samples <- rownames(coldata)[coldata[[group_col]] == group_value]
    sub_counts <- log_counts[, group_samples]
  }
  
  # Calculate pairwise Pearson correlation
  cor_matrix <- cor(sub_counts, method = "pearson")
  
  # Compute average correlation per sample (excluding self-correlation)
  avg_cor <- apply(cor_matrix, 1, function(x) mean(x[names(x) != names(x)[which.max(x)]]))
  
  return(sort(avg_cor))  # Lower = less similar to others
}

#' Compute Euclidean distance to centroid per sample for outlier detection
#'
#' Computes the centroid expression profile across samples and calculates the
#' Euclidean distance of each sample to this centroid, sorted from lowest to highest.
#'
#' @param log_counts A matrix or data.frame of log-scale expression values with samples as columns.
#' @param coldata Sample metadata table with rownames matching sample names.
#' @param group_col Column name used to subset samples when all is FALSE.
#' @param group_value Value in group_col used to select samples when all is FALSE.
#' @param all Logical; if TRUE uses all samples, otherwise subsets by group_col and group_value.
#'
#' @return A named numeric vector of Euclidean distances sorted ascending.
#'
get_euclidean_outliers <- function(log_counts, coldata, group_col, group_value, all = FALSE) {
  if (all) {
    # Use all samples
    sub_counts <- log_counts
  } else {
    # Subset to samples from the specified group
    group_samples <- rownames(coldata)[coldata[[group_col]] == group_value]
    sub_counts <- log_counts[, group_samples]
  }
  
  # Compute the centroid (mean expression per gene)
  centroid <- rowMeans(sub_counts)
  
  # Calculate Euclidean distance of each sample to the centroid
  distances <- apply(sub_counts, 2, function(sample_vec) {
    sqrt(sum((sample_vec - centroid)^2))
  })
  
  return(sort(distances))  # Higher = farther from the group
}
  
#' Plot a sample-to-sample distance heatmap
#'
#' Builds or accepts a DESeqDataSet, applies a chosen transformation, computes a
#' sample-to-sample distance matrix using Euclidean distance or Pearson distance,
#' and visualizes the matrix using pheatmap or an alternative heatmap function.
#'
#' @param counts Either a raw counts matrix or a DESeqDataSet.
#' @param coldata Sample metadata table required when counts is a raw matrix.
#' @param transform Transformation to apply: rlog, vst, or log2.
#' @param distMethod Distance method: euclidean or pearson.
#' @param clusteringMethod Clustering method used for hierarchical clustering.
#' @param usePheatmap Logical; if TRUE uses pheatmap, otherwise uses heatmap.2.
#' @param ... Additional arguments passed through to the plotting function.
#'
#' @return No explicit return value. Produces a heatmap as a side effect.
#'
plotSampleHeatmap <- function(counts,
                              coldata = NULL,
                              transform = c("rlog", "vst", "log2"),
                              distMethod = c("euclidean","pearson"),
                              clusteringMethod = "complete",
                              usePheatmap = TRUE,
                              ...) {
  
  # 1. Build or unwrap the DESeqDataSet
  if (!inherits(counts, "DESeqDataSet")) {
    if (is.null(coldata)) stop("coldata required for raw matrix")
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData   = coldata,
                                  design    = ~ 1)
  } else {
    dds <- counts
  }
  
  # 2. Transform counts
  tr <- match.arg(transform)
  if (tr == "log2") {
    # ensure size factors exist before normalized=TRUE
    if (is.null(sizeFactors(dds)) || all(sizeFactors(dds) == 1)) {
      dds <- estimateSizeFactors(dds)
    }
    mat <- log2(counts(dds, normalized = TRUE) + 1)
  } else if (tr == "rlog") {
    mat <- assay(rlog(dds, blind = TRUE))
  } else {  # vst
    mat <- assay(vst(dds, blind = TRUE))
  }
  
  # 3. Compute distance matrix
  dm <- switch(match.arg(distMethod),
               euclidean = dist(t(mat)),
               pearson   = as.dist(1 - cor(mat, method = "pearson")))
  dmMat <- as.matrix(dm)
  
  # 4. Plot
  pal <- colorRampPalette(rev(brewer.pal(9, "Blues")))(100)
  if (usePheatmap) {
    pheatmap(dmMat,
             clustering_distance_rows = dm,
             clustering_distance_cols = dm,
             clustering_method       = clusteringMethod,
             col                     = pal,
             main                    = "Sample to Sample Distance",
             ...)
  } else {
    heatmap.2(dmMat,
              trace = "none",
              Colv  = as.dendrogram(hclust(dm, method = clusteringMethod)),
              Rowv  = as.dendrogram(hclust(dm, method = clusteringMethod)),
              col   = pal,
              margins = c(10, 10),
              main = "Sample to Sample Distance")
  }
}


#' Plot per-sample expression distributions as violin plots
#'
#' Converts a sample-by-gene expression matrix to long format and plots per-sample
#' distributions using violin plots with an overlaid boxplot. Optional sample order
#' and custom fill palette can be provided.
#'
#' @param expr_mat A matrix or data.frame of expression values with genes as rows and samples as columns.
#' @param sample_order Optional character vector specifying the desired sample order.
#' @param palette Optional vector of fill colors matching the number of samples.
#'
#' @return A ggplot object.
#'
plotSampleViolin <- function(expr_mat, sample_order = NULL, palette = NULL){
  # require packages
  if(!requireNamespace("ggplot2", quietly=TRUE)) stop("Please install ggplot2")
  if(!requireNamespace("tidyr", quietly=TRUE)) stop("Please install tidyr")
  
  # Convert to data.frame
  df <- as.data.frame(expr_mat, stringsAsFactors = FALSE)
  
  # Add gene ID as a column if present
  df$gene <- rownames(df)
  
  # Pivot to long format: columns = sample, expression
  long <- tidyr::pivot_longer(
    df, 
    cols = -gene, 
    names_to  = "Sample", 
    values_to = "Expression"
  )
  
  # If user supplied a sample order, enforce it
  if(!is.null(sample_order)){
    long$Sample <- factor(long$Sample, levels = sample_order)
  }
  
  # Build the violin + boxplot
  p <- ggplot2::ggplot(long, ggplot2::aes_string(x = "Sample", y = "Expression")) +
    ggplot2::geom_violin(trim = TRUE, aes(fill = Sample), show.legend = FALSE) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Sample",
      y = "log2 transformed normalised expression",
      title = "Per-Sample Expression Distributions"
    )
  
  # If user provided a palette of matching length, apply it
  if(!is.null(palette)){
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  
  return(p)
}

#' Plot expression distributions across samples using frequency polygons
#'
#' Creates density-like frequency polygon plots of log-scale expression values
#' per sample, optionally limited to a specified subset of samples.
#'
#' @param log_counts A matrix or data.frame of log2(Expression + 1) values with samples as columns.
#' @param binwidth Bin width used for frequency polygons.
#' @param palette A ggplot2 scale used to color samples.
#' @param y_max Maximum y-axis value shown in the plot.
#' @param subset Optional character vector of sample names to include when all is FALSE.
#' @param all Logical; if TRUE plots all samples, otherwise uses subset.
#'
#' @return A ggplot object.
#'
plotExpressionDistributions <- function(log_counts,
                                        binwidth = 0.2,
                                        palette   = scale_color_viridis_d(option = "turbo"),
                                        y_max     = 3000,
                                        subset    = NULL,
                                        all       = TRUE) {
  # ensure data.frame
  df <- as.data.frame(log_counts)
  
  # determine which samples to include
  samples <- colnames(df)
  if (!all) {
    if (is.null(subset) || length(subset) == 0) {
      stop("When all = FALSE, please provide a non‐empty `subset` of sample names.")
    }
    missing <- setdiff(subset, samples)
    if (length(missing) > 0) {
      stop("These samples were not found: ", paste(missing, collapse = ", "))
    }
    samples_to_plot <- subset
  } else {
    samples_to_plot <- samples
  }
  
  # pivot to long & filter
  df_long <- df %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols      = all_of(samples_to_plot),
      names_to  = "sample",
      values_to = "log2expr"
    )
  
  # build plot
  ggplot(df_long, aes(x = log2expr, color = sample)) +
    geom_freqpoly(binwidth = binwidth, size = 0.8) +
    palette +
    coord_cartesian(ylim = c(0, y_max)) +
    labs(
      title    = "RNA-Seq Expression Distributions by Sample",
      subtitle = if (all) "All samples shown" else paste0("Showing subset: ", paste(samples_to_plot, collapse = ", ")),
      x        = "Log2(Expression + 1)",
      y        = "Number of Genes",
      color    = "Sample"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
}


#' Aggregate samples in a count table by a grouping column
#'
#' Aggregates counts across sets of samples defined by a grouping variable in coldata.
#' Supports summing or averaging counts, and constructs a new coldata table containing
#' only requested metadata columns carried over from the first sample in each group.
#'
#' @param coldata Sample metadata table with rownames matching sample names.
#' @param count_table A raw counts matrix with genes as rows and samples as columns.
#' @param group_col Column name in coldata defining sample groups to aggregate.
#' @param keep_cols Optional character vector of metadata columns to carry over to the new coldata.
#' @param method Aggregation method: sum or average.
#'
#' @return A list with two elements: count_table (aggregated counts) and coldata (aggregated metadata).
#'
aggregate_samples_by <- function(coldata,
                                 count_table,
                                 group_col,
                                 keep_cols = NULL,
                                 method = c("sum","average")) {
  method <- match.arg(method)
  #–– input checks
  if (!group_col %in% colnames(coldata)) {
    stop("'", group_col, "' not found in coldata")
  }
  if (!all(rownames(coldata) %in% colnames(count_table))) {
    stop("Some coldata rows are not columns in count_table")
  }
  if (!is.null(keep_cols) && !all(keep_cols %in% colnames(coldata))) {
    stop("Some keep_cols not in coldata")
  }
  
  groups   <- unique(coldata[[group_col]])
  ngroups  <- length(groups)
  ngenes   <- nrow(count_table)
  
  # prepare outputs
  mat <- matrix(0L, nrow = ngenes, ncol = ngroups,
                dimnames = list(rownames(count_table), NULL))
  
  # Initialize new_meta data frame with ONLY keep_cols (not the group_col)
  if (is.null(keep_cols)) {
    new_meta <- data.frame(matrix(nrow = ngroups, ncol = 0))
  } else {
    new_meta <- data.frame(matrix(nrow = ngroups, ncol = length(keep_cols)))
    colnames(new_meta) <- keep_cols
  }
  
  merged_names <- character(ngroups)
  
  for (i in seq_along(groups)) {
    g       <- groups[i]
    samp    <- rownames(coldata)[coldata[[group_col]] == g]
    
    # carry over metadata - ONLY keep_cols, not group_col
    if (!is.null(keep_cols)) {
      for (col in keep_cols) {
        new_meta[[col]][i] <- as.character(coldata[samp[1], col])
      }
    }
    
    # 1) if only one sample, keep its name
    if (length(samp) == 1) {
      merged <- samp
      counts_agg <- count_table[, samp]
    }
    # 2) if GEO IDs, collapse literally
    else if (all(grepl("^SRR", samp))) {
      merged <- paste(samp, collapse = "_")
      counts_agg <- if (method=="sum") {
        rowSums(count_table[, samp, drop=FALSE])
      } else {
        round(rowMeans(count_table[, samp, drop=FALSE]))
      }
    }
    # 3) else use prefix/ID merging logic
    else {
      prefixes    <- gsub("_S[0-9]+$", "", samp)
      parts       <- strsplit(prefixes, "_")
      common_pref <- parts[[1]][1]
      ids         <- sapply(parts, function(p) if (length(p)>1) p[2] else "")
      merged      <- paste0(common_pref, "_", paste(ids, collapse = "_"))
      counts_agg  <- if (method=="sum") {
        rowSums(count_table[, samp, drop=FALSE])
      } else {
        round(rowMeans(count_table[, samp, drop=FALSE]))
      }
    }
    
    # fill outputs
    merged_names[i]     <- merged
    mat[, i]            <- counts_agg
  }
  
  # finalize
  colnames(mat)         <- merged_names
  rownames(new_meta)    <- merged_names
  new_counts_df         <- as.data.frame(mat, stringsAsFactors = FALSE)
  
  list(count_table = new_counts_df, coldata = new_meta)
}

# Convenience wrappers:
sum_samples_by <- function(coldata, count_table, group_col, keep_cols = NULL) {
  aggregate_samples_by(coldata, count_table, group_col, keep_cols, method = "sum")
}

average_samples_by <- function(coldata, count_table, group_col, keep_cols = NULL) {
  aggregate_samples_by(coldata, count_table, group_col, keep_cols, method = "average")
}


