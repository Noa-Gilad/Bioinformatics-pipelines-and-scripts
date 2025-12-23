# ==============================================================================
# Author: Noa Gilad
# Title: General DESeq2 Pipeline Script
# Description: 
#     This R script implements a general-purpose DESeq2 differential expression 
#     analysis pipeline. It assumes input counts and metadata tables and provides 
#     an interface to generate PCA plots, differential gene lists, and other 
#     common outputs.
#
# Note:
#     This script is generic and not tied to a specific project. To use it, 
#     provide a raw counts matrix and matching metadata with appropriate conditions.
#
# Dependencies: DESeq2, ggplot2, pheatmap, dplyr, etc.
# ==============================================================================

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Load required libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# tidy-verse style data wrangling ------------------------------------------------
library(dplyr)            # %>%, mutate(), arrange(), join helpers, piping
library(readr)            # read_csv(), write_tsv(), fast Delimiter-based I/O
library(tibble)           # tibble(), as_tibble(); keeps data frames tidy

# visualisation helpers ---------------------------------------------------------
library(ggplot2)          # grammar of graphics base layer
library(ggpubr)           # publication-ready wrappers (stat_compare_means, etc.)
library(RColorBrewer)     # qualitative & sequential colour palettes
library(ggforce)          # geom_mark_ellipse() used inside custom do_pca()
library(gridExtra)        # arrangeGrob(), grid.arrange() for multi-plot figures

# core RNA-seq / batch-effect workflow -----------------------------------------
library(DESeq2)           # DESeqDataSet, DESeq(), counts(), results(), vst()
library(sva)              # svaseq() for surrogate-variable analysis (SVA)
library(limma)            # removeBatchEffect() for batch / SVA correction

# enrichment / annotation -------------------------------------------------------
library(clusterProfiler)  # enrichGO(), enrichKEGG(), dotplot() for GO/KEGG
library(org.Hs.eg.db)     # human gene→GO mapping (OrgDb) used by clusterProfiler
library(AnnotationDbi)

# misc utilities ----------------------------------------------------------------
library(matrixStats)      # rowVars(), colMedians()—fast matrix statistics
library(reshape2)         # melt() for reshaping matrices to long format
library(glue)             # glue() string interpolation for filenames, messages



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' Prepare and run a DESeq2 analysis object
#'
#' Creates a DESeqDataSet from a raw counts matrix and sample metadata, relevels
#' the specified condition factor to a reference level, then runs the DESeq
#' model fitting procedure.
#'
#' @param counts A raw counts matrix (genes by samples).
#' @param coldata A sample metadata table with rownames matching sample names.
#' @param refer Reference level for the condition factor.
#' @param des Design formula for DESeq2 (for example, ~ condition).
#' @param condition Column name in coldata indicating the experimental condition.
#'
#' @return A DESeqDataSet after running DESeq.
#'
# load data
prepare_data <- function(counts, coldata, refer, des, condition) {
  coldata[[condition]] <- relevel(coldata[[condition]], ref = refer)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = des)
  dds <- DESeq(dds)
  return(dds)
} 

#' Extract DESeq2 results for a specific contrast
#'
#' Computes differential expression results for a requested contrast, sorts by
#' adjusted p-value, optionally writes the result table to CSV, and returns the
#' results object.
#'
#' @param dds A DESeqDataSet that has already been processed with DESeq.
#' @param toCompare The numerator level in the contrast (for example, treated).
#' @param baseLevel The denominator/base level in the contrast (for example, control).
#' @param lfc Log2 fold-change threshold used by DESeq2 for thresholded testing.
#' @param alphaLim Significance level used by DESeq2 for multiple testing.
#' @param condition Column name in colData(dds) used as the contrast factor.
#' @param writeFile Logical; if TRUE, writes "<toCompare>_vs_<baseLevel>_DGE.csv".
#'
#' @return A DESeqResults object sorted by adjusted p-value.
#'
# look at the results of dds
results_data <- function(dds, toCompare, baseLevel, lfc, alphaLim, condition, writeFile) {
  res <- results(dds, lfcThreshold = lfc, alpha = alphaLim, contrast = c(condition, toCompare, baseLevel)) 
  head(results(dds, tidy=TRUE))
  res <- res[order(res$padj),]
  head(res,10)
  summary(res)
  if (writeFile == TRUE) {
    write.csv(as.data.frame(res), file=glue('{toCompare}_vs_{baseLevel}_DGE.csv'))
  }
  return(res)
}

#' PCA on top-variable genes after low-count filtering and VST
#'
#' Performs low-count filtering, variance-stabilizing transformation, optional
#' batch and/or surrogate variable removal, selects top variable genes, runs PCA,
#' and returns either the PCA data frame or a ggplot.
#'
#' @param dds A DESeqDataSet.
#' @param first_pc Principal component to plot on the x-axis.
#' @param sec_pc Principal component to plot on the y-axis.
#' @param intgroup Column name(s) in colData used for grouping/coloring.
#' @param ntop Number of most-variable genes used for PCA.
#' @param min_count Minimum raw count to consider a gene expressed in a sample.
#' @param min_samples Minimum number of samples meeting min_count (defaults to smallest group size).
#' @param returnData Logical; if TRUE returns the PCA data frame instead of the plot.
#' @param pca_colors Colors used for group mapping in the plot.
#' @param batch Logical; if TRUE applies batch removal.
#' @param batch_rm Batch variable vector used for batch removal.
#' @param sva Logical; if TRUE applies surrogate variable removal.
#' @param sva_rm Matrix of surrogate variables aligned to samples.
#' @param design Design matrix required for batch/SVA removal.
#'
#' @return Either a ggplot PCA plot, or a data.frame if returnData is TRUE.
#'
# PCA on top-variable genes with low-count filtering (pre-VST)
do_pca <- function(dds,
                   first_pc    = 1,
                   sec_pc      = 2,
                   intgroup,
                   ntop        = 500,
                   min_count   = 10,
                   min_samples = NULL,
                   returnData  = FALSE,
                   pca_colors  = c("black","blue"),
                   batch       = FALSE,
                   batch_rm    = NULL,
                   sva         = FALSE,
                   sva_rm      = NULL,    # your svseq$sv matrix
                   design      = NULL) {
  
  # 1) default min_samples
  if (is.null(min_samples)) {
    min_samples <- min(table(colData(dds)[, intgroup[1]]))
  }
  
  # 2) low-count filter
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  if (sum(keep) < ntop) {
    warning("Only ", sum(keep), " genes pass. Reducing ntop to ", sum(keep), ".")
    ntop <- sum(keep)
  }
  dds <- dds[keep, ]
  
  # 3) VST
  vsd <- vst(dds)
  
  # 4) batch/SVA removal
  if (batch || sva) {
    if (batch && is.null(batch_rm)) {
      stop("When batch=TRUE, supply batch_rm.")
    }
    if (sva && is.null(sva_rm)) {
      stop("When sva=TRUE, supply sva_rm (svseq$sv).")
    }
    if (is.null(design)) {
      stop("When batch or sva=TRUE, supply design (model.matrix).")
    }
    
    # reorder your sva matrix if it has rownames
    if (sva && !is.null(rownames(sva_rm))) {
      sva_rm <- sva_rm[colnames(vsd), , drop = FALSE]
    }
    
    assay(vsd) <- removeBatchEffect(
      x          = assay(vsd),
      batch      = if (batch) batch_rm else NULL,
      covariates = if (sva) sva_rm else NULL,
      design     = design
    )
  }
  
  # 5) pick top-variable genes
  mat <- assay(vsd)
  rv  <- rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(ntop)]
  
  # 6) PCA
  pca        <- prcomp(t(mat[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  
  # 7) grouping
  meta  <- as.data.frame(colData(vsd)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(meta, 1, paste, collapse = ":"))
  } else {
    colData(vsd)[, intgroup]
  }
  
  # 8) assemble
  df <- data.frame(
    PC1   = pca$x[, first_pc],
    PC2   = pca$x[, sec_pc],
    group = group,
    meta,
    name  = colnames(vsd)
  )
  attr(df, "percentVar") <- percentVar[c(first_pc, sec_pc)]
  if (returnData) return(df)
  
  # 9) plot
  ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    geom_mark_ellipse(aes(fill = group), alpha = 0.2) +
    geom_text(aes(label = name), vjust = -0.5, hjust = 1) +
    scale_color_manual(values = pca_colors) +
    scale_fill_manual(values  = pca_colors) +
    xlab(sprintf("PC%d (%.1f%%)", first_pc,  100 * percentVar[first_pc])) +
    ylab(sprintf("PC%d (%.1f%%)", sec_pc,    100 * percentVar[sec_pc])) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = 0.20)) +
    theme_minimal(base_size = 12)
}

#' Plot per-sample VST expression distributions as boxplots
#'
#' Converts a VST-transformed object to long format and returns a ggplot boxplot
#' of expression values per sample.
#'
#' @param vst_object A VST-transformed object (output of vst).
#'
#' @return A ggplot object.
#'
plot_vst_boxplot <- function(vst_object) {
  # Extract the expression matrix from the VST object
  df <- as.data.frame(assay(vst_object))
  
  # Add gene names as a column (the row names of the matrix)
  df$Gene <- rownames(df)
  
  # Melt the data frame into long format using reshape2::melt
  df_long <- melt(df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
  
  # Set the factor levels for Sample according to the original column order (excluding the Gene column)
  df_long$Sample <- factor(df_long$Sample, levels = colnames(df)[-ncol(df)])
  
  # Create the box plot without log-transformation because VST is already on a log scale
  p <- ggplot(df_long, aes(x = Sample, y = Expression)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = "VST Expression Distribution Across Samples",
         x = "Samples",
         y = "VST Expression Values") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(p)
}


#' Extract significantly differentially expressed genes as an ID table
#'
#' Filters a results table for up- or down-regulated genes using log2 fold-change
#' and adjusted p-value thresholds, returning only an ID column.
#'
#' @param res A results object or data.frame containing log2FoldChange and padj.
#' @param direction "up" for upregulated, otherwise treated as downregulated.
#' @param lfc Log2 fold-change threshold used for filtering.
#' @param alphaLim Adjusted p-value cutoff.
#'
#' @return A data.frame with a single column named ID containing passing gene rownames.
#'
# Extracts differentialy expressed genes (up or down)
extract_DEG <- function(res, direction = "up", lfc, alphaLim) {
  res <- na.omit(res)  # Remove rows with missing values (NA)
  
  # Filter for upregulated or downregulated genes based on the direction
  if (direction == "up") {
    genes <- res[res$log2FoldChange > lfc & res$padj < alphaLim, ]
  } else {
    genes <- res[res$log2FoldChange < -lfc & res$padj < alphaLim, ]
  }
  
  return(data.frame(ID = rownames(genes)))
}

#' Extract DE genes and return the full filtered results table
#'
#' Filters a results table for up- or down-regulated genes using log2 fold-change
#' and adjusted p-value thresholds, sorts by adjusted p-value, adds gene symbols,
#' optionally writes the filtered table to CSV, and returns it.
#'
#' @param res A results object or data.frame containing log2FoldChange and padj.
#' @param direction "up" for upregulated, otherwise treated as downregulated.
#' @param lfc Log2 fold-change threshold used for filtering.
#' @param alphaLim Adjusted p-value cutoff.
#' @param writeFile Logical; if TRUE writes "<title>.csv".
#' @param title Output filename prefix used when writeFile is TRUE.
#'
#' @return A data.frame of filtered genes sorted by adjusted p-value with an added SYMBOL column.
#'
# Extracts differentialy expressed genes (up or down)
extract_DEG_deseq2_data <- function(res, direction = "up", lfc, alphaLim, writeFile, title) {
  res <- na.omit(res)  # Remove rows with missing values (NA)
  
  # Filter for upregulated or downregulated genes based on the direction
  if (direction == "up") {
    genes <- res[res$log2FoldChange > lfc & res$padj < alphaLim, ]
  } else {
    genes <- res[res$log2FoldChange < -lfc & res$padj < alphaLim, ]
  }
  
  # Sort by padj in ascending order (most significant first)
  genes <- genes[order(genes$padj), ]
  
  genes$SYMBOL <- mapIds(org.Hs.eg.db, keys = rownames(genes), keytype = "ENSEMBL", column = "SYMBOL")
  
  if (writeFile == TRUE) {
    write.csv(as.data.frame(genes), file=glue('{title}.csv'))
  }
  
  return(genes)
}


#' Perform Gene Ontology enrichment analysis and generate a barplot
#'
#' Runs GO enrichment on a set of gene IDs, handles errors gracefully, and returns
#' both a results table and a barplot of top categories.
#'
#' @param gene_ids Vector of gene identifiers to test for enrichment.
#' @param org_db Organism annotation database.
#' @param id_type Key type for gene_ids.
#' @param ont GO ontology: BP, MF, or CC.
#' @param universe Optional background gene list.
#' @param show_category Number of categories to show in the barplot.
#' @param p_cutoff P-value cutoff for enrichment.
#' @param q_cutoff FDR cutoff for enrichment.
#'
#' @return A list with results_df (data.frame) and plot (barplot object or NULL).
#'
perform_GO_enrichment <- function(gene_ids, 
                                  org_db = org.Hs.eg.db,  # Default to human database
                                  id_type = "ENSEMBL",    # Default to ENSEMBL IDs
                                  ont = "BP",             # Default to Biological Process
                                  universe = NULL,
                                  show_category = 15,     # Number of categories to show 
                                  p_cutoff = 0.05,        # p-value cutoff
                                  q_cutoff = 0.2) {       # q-value cutoff
  
  # Perform GO enrichment analysis
  go_results <- tryCatch({
    enrichGO(
      gene         = gene_ids,     # List of gene IDs to analyze
      OrgDb        = org_db,       # Organism database for GO terms
      keyType      = id_type,      # Input gene ID type
      ont          = ont,          # Ontology category: "BP", "MF", or "CC"
      universe     = universe,
      pAdjustMethod = "BH",        # Method for adjusting p-values (Benjamini-Hochberg)
      pvalueCutoff = p_cutoff,     # P-value cutoff for significant terms
      qvalueCutoff = q_cutoff      # q-value cutoff (FDR-adjusted p-value)
    )
  }, error = function(e) {
    message("Error in GO enrichment: ", e$message)
    return(NULL)
  })
  
  # Check if enrichment was successful
  if (is.null(go_results) || nrow(go_results) == 0) {
    message("No significant GO terms found.")
    return(list(
      results_df = data.frame(),
      plot = NULL
    ))
  }
  
  # Convert the GO enrichment results to a data frame
  go_results_df <- as.data.frame(go_results)
  
  # Create a bar plot for the top GO categories
  go_plot <- barplot(go_results, showCategory = min(show_category, nrow(go_results_df)))
  
  # Return a list containing the results data frame and the bar plot
  return(list(
    results_df = go_results_df,  # Data frame of GO enrichment results
    plot = go_plot               # Bar plot of the top GO categories
  ))
}

#' Create a customized MA plot from results
#'
#' Builds an MA plot using mean expression versus log2 fold change, colors points
#' by significance using FDR and fold-change thresholds, and optionally highlights
#' a provided gene list.
#'
#' @param res A results data.frame with baseMean, log2FoldChange, and padj columns.
#' @param gene_list Optional data.frame of genes to highlight; must contain an ID column.
#' @param highlight_color Color used for highlighted genes.
#' @param up_color Color for significant up-regulated genes.
#' @param down_color Color for significant down-regulated genes.
#' @param ns_color Color for non-significant genes.
#' @param fdr Adjusted p-value cutoff.
#' @param fc Fold-change threshold on the linear scale.
#' @param size Point size for the base scatter layer.
#' @param alpha Point alpha for the base scatter layer.
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#'
#' @return A ggplot MA plot.
#'
do_MA_plot <- function(res, gene_list = NULL,  # Optional gene list to color 
                       highlight_color = "purple", # Color for highlighted genes
                       up_color = "#B31B21", # Red for up-regulated
                       down_color = "#007FFF", # Blue for down-regulated
                       ns_color = "darkgray", # Gray for non-significant
                       fdr = 0.05, 
                       fc = 1.5, 
                       size = 1, 
                       alpha = 1,
                       main = "MA Plot",
                       xlab = "Log2 mean expression", 
                       ylab = "Log2 fold change") {
  
  # Check data format
  req_cols <- c("baseMean", "log2FoldChange", "padj")
  if(!all(req_cols %in% colnames(res))) {
    stop("Input data must contain: ", paste(req_cols, collapse = ", "))
  }
  
  # Create data frame for plotting
  data <- data.frame(
    name = rownames(res),
    mean = log2(res$baseMean + 1),  # Log2 transform if not already
    lfc = res$log2FoldChange,
    padj = res$padj
  )
  
  # Determine significance status
  sig <- rep("NS", nrow(data))  # Default is non-significant (NS)
  sig[which(data$padj <= fdr & data$lfc < 0 & abs(data$lfc) >= log2(fc))] <- "Down"  # Down-regulated
  sig[which(data$padj <= fdr & data$lfc > 0 & abs(data$lfc) >= log2(fc))] <- "Up"  # Up-regulated
  data$sig <- factor(sig, levels = c("Up", "Down", "NS"))
  
  # Add count labels to the legend
  data$sig_label <- factor(
    sig,
    levels = c("Up", "Down", "NS"),
    labels = c(
      paste0("Up: ", sum(sig == "Up")),
      paste0("Down: ", sum(sig == "Down")),
      "NS"
    )
  )
  
  # Create base plot
  p <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(color = sig_label), size = size, alpha = alpha) +
    scale_color_manual(values = c(up_color, down_color, ns_color)) +
    labs(x = xlab, y = ylab, title = main, color = "") +
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), 
               linetype = c(1, 2, 2),
               color = c("black", "black", "black")) +
    theme_minimal()
  
  # If gene list is provided, highlight those genes
  if (!is.null(gene_list) && nrow(gene_list) > 0) {
    # Check if gene_list has the right format
    if (!("ID" %in% colnames(gene_list))) {
      stop("gene_list must have a column named 'ID'")
    }
    
    # Find genes in the list that are in our data
    highlight_genes <- subset(data, name %in% gene_list$ID)
    
    if (nrow(highlight_genes) > 0) {
      p <- p + 
        geom_point(data = highlight_genes, color = highlight_color, size = 2) 
    }
  }
  
  return(p)
}


#' Apply SVA correction in a DESeq2 workflow
#'
#' Identifies surrogate variables using normalized counts and user-provided full
#' and null models, adds the surrogate variables to the dataset, updates the design,
#' runs DESeq again, and returns the corrected object and surrogate variable values.
#'
#' @param dds DESeq2 object with the raw data.
#' @param full_model Full model formula.
#' @param null_model Null model formula.
#' @param n_sv Number of surrogate variables to use.
#' @param low_count_filter Threshold for filtering low-expressed genes.
#'
#' @return A list containing the corrected DESeq2 object and the surrogate variables.
#'
apply_sva_correction <- function(dds,
                                 full_model,
                                 null_model,
                                 n_sv,
                                 low_count_filter = 1) {
  
  # Extract normalized counts data
  dat <- counts(dds, normalized = TRUE)
  
  # Filter low-expressed genes
  idx <- rowMeans(dat) > low_count_filter
  dat <- dat[idx, ]
  
  # Create model matrices
  mod <- model.matrix(full_model, colData(dds))
  mod0 <- model.matrix(null_model, colData(dds))
  
  # --- SVA Analysis ---
  cat("Performing SVA with n.sv =", n_sv, "\n")
  svseq <- svaseq(dat, mod, mod0, n.sv = n_sv)
  
  # Create a new DESeq2 object and add the SVs
  ddssva <- dds
  for (i in 1:n_sv) {
    ddssva[[paste0("SV", i)]] <- svseq$sv[, i]
  }
  
  # Get the original design formula as a character string
  original_design <- as.character(full_model)[2]  # Remove the ~ part
  
  # Update the design to include SVs
  sv_terms <- paste0("SV", 1:n_sv, collapse = " + ")
  new_design <- formula(paste("~", sv_terms, "+", original_design))
  design(ddssva) <- new_design
  
  # Run DESeq on the SVA-corrected object
  ddssva <- DESeq(ddssva)
  
  # Create result list
  result <- list(
    ddssva = ddssva,     # The DESeq object with SVA included in the design
    sv_values = svseq$sv # The actual surrogate variables
  )
  
  # Return the results
  return(result)
}

#' Generate PCA plots across a range of surrogate variable counts
#'
#' Iteratively runs SVA correction for a range of surrogate variable counts, then
#' runs PCA on each corrected dataset with optional batch removal and optional plot saving.
#'
#' @param dds A DESeqDataSet.
#' @param full_model Full model formula for SVA.
#' @param null_model Null model formula for SVA.
#' @param vars_to_preserve Variables to preserve in the design matrix used for correction.
#' @param min_nsv Minimum number of surrogate variables to test.
#' @param max_nsv Maximum number of surrogate variables to test.
#' @param batch_var Optional column name in colData treated as an explicit batch factor.
#' @param plots_per_grid Number of PCA plots per saved grid image.
#' @param save_plots Logical; if TRUE saves plot grids as PNG files.
#' @param output_prefix Prefix for saved PNG filenames.
#' @param plot_width Width for saved grid images.
#' @param plot_height Height for saved grid images.
#' @param pca_colors Colors passed to do_pca.
#' @param pca_intgroup Column(s) used for PCA grouping; defaults to first variable in full_model.
#' @param low_count_filter Low-count filter passed to apply_sva_correction.
#'
#' @return A list of PCA ggplot objects, one per tested surrogate variable count.
#'
visualize_nsv_choice <- function(dds, 
                                 full_model,
                                 null_model,
                                 vars_to_preserve,
                                 min_nsv = 1, 
                                 max_nsv = 30,
                                 batch_var = NULL,
                                 plots_per_grid = 5,
                                 save_plots = TRUE,
                                 output_prefix = "PCA_plots_nsv",
                                 plot_width = 20,
                                 plot_height = 5,
                                 pca_colors = c("black", "blue"),
                                 pca_intgroup = NULL,
                                 low_count_filter = 1) {
  
  # --- sanity checks ---
  if (max_nsv < min_nsv) stop("max_nsv must be >= min_nsv")
  if (is.null(pca_intgroup)) {
    vars_in_model <- all.vars(full_model)
    if (length(vars_in_model) < 1) stop("No pca_intgroup and full_model has no variables")
    pca_intgroup <- vars_in_model[1]
  }
  
  # prepare formula for preservation
  design_formula <- as.formula(
    paste("~", paste(vars_to_preserve, collapse = " + "))
  )
  
  pca_plots <- vector("list", length = max_nsv - min_nsv + 1)
  idx        <- 1
  
  for (n.sv in seq(min_nsv, max_nsv)) {
    message("Processing n.sv = ", n.sv)
    
    # 1) get SVA correction
    sva_out <- apply_sva_correction(
      dds               = dds,
      full_model        = full_model,
      null_model        = null_model,
      n_sv              = n.sv,
      low_count_filter  = low_count_filter
    )
    ddssva  <- sva_out$ddssva
    sv_vals <- sva_out$sv_values  # matrix (samples × n.sv)
    
    # 2) build design matrix to preserve vars_to_preserve
    design_mat <- model.matrix(design_formula, colData(ddssva))
    
    # 3) extract batch factor if requested
    if (!is.null(batch_var)) {
      if (! batch_var %in% colnames(colData(ddssva))) {
        stop("batch_var '", batch_var, "' not found in colData(dds)")
      }
      batch_flag <- TRUE
      batch_rm   <- colData(ddssva)[[batch_var]]
    } else {
      batch_flag <- FALSE
      batch_rm   <- NULL
    }
    
    # 4) call do_pca()
    p <- do_pca(
      ddssva,
      intgroup   = pca_intgroup,
      batch      = batch_flag,
      batch_rm   = batch_rm,
      sva        = TRUE,
      sva_rm     = sv_vals,
      design     = design_mat,
      pca_colors = pca_colors
    ) + ggtitle(paste("PCA: n.sv =", n.sv))
    
    pca_plots[[idx]] <- p
    idx <- idx + 1
  }
  
  # --- save grids of plots ---
  if (save_plots) {
    num_plots <- length(pca_plots)
    num_grids <- ceiling(num_plots / plots_per_grid)
    
    for (g in seq_len(num_grids)) {
      start_i <- (g - 1) * plots_per_grid + 1
      end_i   <- min(g * plots_per_grid, num_plots)
      grid    <- gridExtra::arrangeGrob(
        grobs = pca_plots[start_i:end_i],
        ncol  = min(plots_per_grid, end_i - start_i + 1),
        top   = paste0("PCA: n.sv ", start_i + min_nsv - 1,
                       " to ",    end_i   + min_nsv - 1)
      )
      fname <- sprintf("%s_%d_to_%d.png",
                       output_prefix,
                       start_i + min_nsv - 1,
                       end_i   + min_nsv - 1)
      ggsave(filename = fname, plot = grid,
             width = plot_width, height = plot_height)
    }
  }
  
  return(pca_plots)
}


