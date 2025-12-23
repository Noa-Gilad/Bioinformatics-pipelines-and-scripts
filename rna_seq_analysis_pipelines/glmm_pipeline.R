# ==============================================================================
# Author: Noa Gilad
# Title: GLMMSeq Differential Expression Pipeline
# Description:
#   Utility functions for differential expression analysis using generalized
#   linear mixed models with glmmSeq. Includes dispersion and size factor
#   estimation, result formatting for two-level factors, and extraction of
#   significantly differentially expressed genes.
#
# Notes:
#   - Designed for count-based RNA-seq data with fixed effects and a random factor.
#   - Functions are project-agnostic and assume ENSEMBL gene identifiers.
#
# Dependencies: DESeq2, glmmSeq, limma, dplyr, ggplot2, org.Hs.eg.db
# ==============================================================================

#  ~~~~~~~~~~~~~~~~~~~~~~~ Load necessary libraries ~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(readr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(glue)
library(glmmSeq)
library(limma)
library(org.Hs.eg.db)
set.seed(1234)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Estimate gene-wise dispersions using DESeq2
#'
#' Estimates dispersion values for each gene by fitting a DESeq2 model with
#' an intercept-only design. These dispersions can be reused in downstream
#' generalized linear mixed model analyses.
#'
#' @param counts A raw counts matrix with genes as rows and samples as columns.
#' @param metadata Sample metadata table with rownames matching sample names.
#'
#' @return A named numeric vector of gene-wise dispersion estimates.
#'
estimate_disp <- function(counts, metadata) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ 1)
  dds <- DESeq(dds)
  dispersions <- setNames(dispersions(dds), rownames(counts))
  return(dispersions)
}


#' Estimate size factors from a counts matrix
#'
#' Computes normalization size factors directly from a raw counts matrix.
#' This is typically used to normalize for sequencing depth prior to
#' generalized linear mixed model fitting.
#'
#' @param counts A raw counts matrix with genes as rows and samples as columns.
#'
#' @return A numeric vector of size factors, one per sample.
#'
estimate_size_factor <- function(counts) {
  return(estimateSizeFactorsForMatrix(counts))
}

#' Format glmmSeq results for a two-level factor
#'
#' Extracts coefficients, standard errors, mean expression values, and
#' statistical significance metrics from a glmmSeq object for a specified
#' two-level categorical variable.
#'
#' @param object A fitted glmmSeq model object.
#' @param var Metadata column name corresponding to the tested variable.
#' @param target Target level of the variable compared to the reference level.
#'
#' @return A data.frame containing effect sizes, standard errors, mean expression,
#'         p-values, and q-values for each gene.
#'
formatGlmmSeqResults <- function(object, var, target) {
  # Metadata and predictions
  md    <- object@modelData
  preds <- object@predict
  
  # Log2 fold-change with pseudocount = 1
  betaCond <- object@stats$coef[,paste0(var, target)]
  
  # Extract statistics
  # SE is stored per-coefficient, named e.g. "ageO"
  se_vals <- object@stats$stdErr[, paste0(var, target)]
  # p- and q-values are stored per-variable, named e.g. "age"
  p_vals  <- object@stats$pvals[, var]
  q_vals  <- object@stats$qvals[, var]
  meanExp  <- object@stats$res[ , "meanExp"]
  
  # Assemble results
  results <- data.frame(
    gene            = rownames(preds),
    coef            = betaCond,
    SE              = se_vals,
    meanExp         = meanExp,
    p_value         = p_vals,
    q_value         = q_vals,
    stringsAsFactors = FALSE
  )
  rownames(results) <- results$gene
  results$gene <- NULL
  return(results)
}

#' Extract differentially expressed genes as an ID list
#'
#' Filters glmmSeq results to identify significantly up- or down-regulated
#' genes based on an effect size threshold and q-value cutoff.
#'
#' @param res A results data.frame produced by formatGlmmSeqResults.
#' @param direction "up" for positive effects or any other value for negative effects.
#' @param beta Effect size threshold used for filtering.
#' @param alphaLim q-value cutoff for statistical significance.
#'
#' @return A data.frame with a single column named ID containing gene identifiers.
#'
extract_DE_genes <- function(res, direction = "up", beta, alphaLim) {
  res <- na.omit(res)  # Remove rows with missing values (NA)
  
  # Filter for upregulated or downregulated genes based on the direction
  if (direction == "up") {
    genes <- res[res$coef > beta & res$q_value < alphaLim, ]
  } else {
    genes <- res[res$coef < -beta & res$q_value < alphaLim, ]
  }
  genes <- genes[order(genes$q_value), ]
  return(data.frame(ID = rownames(genes)))
}

#' Extract differential expression results with annotations
#'
#' Filters glmmSeq results for significantly differentially expressed genes,
#' adds gene symbol annotations, optionally writes the results to a CSV file,
#' and returns the filtered table.
#'
#' @param res A results data.frame produced by formatGlmmSeqResults.
#' @param direction "up" for positive effects or any other value for negative effects.
#' @param beta Effect size threshold used for filtering.
#' @param alphaLim q-value cutoff for statistical significance.
#' @param writeFile Logical; if TRUE writes the results table to disk.
#' @param title Output filename prefix used when writeFile is TRUE.
#'
#' @return A data.frame of filtered genes sorted by q-value with an added SYMBOL column.
#'
extract_DE_data <- function(res, direction = "up", beta, alphaLim, writeFile, title) {
  res <- na.omit(res)  # Remove rows with missing values (NA)
  
  # Filter for upregulated or downregulated genes based on the direction
  if (direction == "up") {
    genes <- res[res$coef > beta & res$q_value < alphaLim, ]
  } else {
    genes <- res[res$coef < -beta & res$q_value < alphaLim, ]
  }
  
  # Sort by padj in ascending order (most significant first)
  genes <- genes[order(genes$q_value), ]
  genes$SYMBOL <- mapIds(org.Hs.eg.db, keys = rownames(genes), keytype = "ENSEMBL", column = "SYMBOL")
  
  if (writeFile == TRUE) {
    write.csv(as.data.frame(genes), file=glue('{title}.csv'))
  }
  
  return(genes)
}


