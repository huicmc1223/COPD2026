###############################################################################
# batch_correction.R — Unified batch correction functions
# Usage: source(here::here("scripts/utils/batch_correction.R"))
###############################################################################

suppressPackageStartupMessages({
  library(sva)
  library(limma)
})

#' Batch correction for raw counts (DESeq2-compatible)
#'
#' Uses ComBat-seq to remove batch effects while preserving count nature.
#' @param counts Integer count matrix (genes × samples).
#' @param batch Factor vector of batch labels (same length as ncol(counts)).
#' @param group Factor vector of biological group (condition).
#' @return Adjusted count matrix (integer).
combat_seq_correct <- function(counts, batch, group = NULL) {
  message("── ComBat-seq batch correction ──")
  message("  Batches: ", paste(levels(factor(batch)), collapse = ", "))
  message("  Samples: ", ncol(counts))
  
  adjusted <- ComBat_seq(counts  = as.matrix(counts),
                         batch   = batch,
                         group   = group,
                         full_mod = TRUE)
  
  message("  ✔ ComBat-seq correction completed")
  return(adjusted)
}

#' Batch correction for normalised expression (microarray / vst)
#'
#' Uses ComBat to remove batch effects from log-transformed expression.
#' @param expr Numeric matrix (genes × samples), log2-scale.
#' @param batch Factor vector of batch labels.
#' @param mod Model matrix for biological variables to preserve (optional).
#' @return Adjusted expression matrix.
combat_correct <- function(expr, batch, mod = NULL) {
  message("── ComBat batch correction ──")
  message("  Batches: ", paste(levels(factor(batch)), collapse = ", "))
  
  adjusted <- ComBat(dat   = as.matrix(expr),
                     batch = batch,
                     mod   = mod,
                     par.prior = TRUE,
                     prior.plots = FALSE)
  
  message("  ✔ ComBat correction completed")
  return(adjusted)
}

#' Remove batch effect using limma (for visualisation / downstream, not DEG)
#'
#' @param expr Numeric matrix (genes × samples), log2-scale.
#' @param batch Factor vector of batch labels.
#' @param design Model matrix for biological effects to preserve.
#' @return Adjusted expression matrix.
limma_remove_batch <- function(expr, batch, design = NULL) {
  message("── limma::removeBatchEffect ──")
  adjusted <- removeBatchEffect(as.matrix(expr), batch = batch, design = design)
  message("  ✔ limma batch correction completed")
  return(adjusted)
}

#' Harmony integration for Seurat objects (scRNA-seq)
#'
#' @param seu Merged Seurat object with metadata column for batch.
#' @param batch_var Name of metadata column with batch labels (default "dataset").
#' @param dims Number of PCA dimensions to use (default 30).
#' @param theta Diversity clustering penalty (default 2).
#' @return Seurat object with Harmony-corrected embeddings.
harmony_integrate <- function(seu, batch_var = "dataset", dims = 30, theta = 2) {
  if (!requireNamespace("harmony", quietly = TRUE)) stop("harmony not installed")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat not installed")
  
  message("── Harmony integration ──")
  message("  Batch variable: ", batch_var)
  message("  Levels: ", paste(unique(seu@meta.data[[batch_var]]), collapse = ", "))
  
  seu <- harmony::RunHarmony(seu,
                             group.by.vars = batch_var,
                             dims.use      = 1:dims,
                             theta         = theta,
                             plot_convergence = FALSE)
  
  # Re-run UMAP on corrected embeddings
  seu <- Seurat::RunUMAP(seu, reduction = "harmony", dims = 1:dims)
  seu <- Seurat::FindNeighbors(seu, reduction = "harmony", dims = 1:dims)
  seu <- Seurat::FindClusters(seu, resolution = 0.5)
  
  message("  ✔ Harmony integration completed — ",
          length(unique(Seurat::Idents(seu))), " clusters")
  return(seu)
}

#' Diagnostic PCA plot: before vs after batch correction
#'
#' @param expr_before Expression matrix before correction.
#' @param expr_after Expression matrix after correction.
#' @param batch Factor vector of batch labels.
#' @param condition Factor vector of conditions (biological).
#' @param title_prefix Prefix for plot title.
#' @return patchwork combined plot.
plot_batch_pca <- function(expr_before, expr_after, batch, condition,
                           title_prefix = "") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 needed")
  
  .pca_df <- function(expr, label) {
    pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)
    pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    data.frame(
      PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      batch = batch, condition = condition,
      stage = label,
      pc1_label = paste0("PC1 (", pct[1], "%)"),
      pc2_label = paste0("PC2 (", pct[2], "%)")
    )
  }
  
  df <- rbind(.pca_df(expr_before, "Before"), .pca_df(expr_after, "After"))
  df$stage <- factor(df$stage, levels = c("Before", "After"))
  
  p_batch <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = batch)) +
    ggplot2::geom_point(size = 1, alpha = 0.7) +
    ggplot2::facet_wrap(~stage, scales = "free") +
    ggplot2::labs(title = paste0(title_prefix, " — coloured by Batch")) +
    theme_nature()
  
  p_cond <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = condition)) +
    ggplot2::geom_point(size = 1, alpha = 0.7) +
    ggplot2::facet_wrap(~stage, scales = "free") +
    ggplot2::labs(title = paste0(title_prefix, " — coloured by Condition")) +
    theme_nature()
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(p_batch / p_cond)
  }
  return(list(batch = p_batch, condition = p_cond))
}

message("✔ batch_correction.R loaded")
