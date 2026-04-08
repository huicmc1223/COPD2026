###############################################################################
# copd_marker_interaction.R — TMEM164 interaction with known COPD markers
#
# Step 13: Correlation + co-expression + WGCNA module analysis
#
# Markers: MUC5B, MMP7, SPP1, HHIP, CC16 (SCGB1A1)
#
# Run: Rscript scripts/03_copd_analysis/copd_marker_interaction.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(pheatmap)
})

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bulk-level correlation (across all COPD bulk/array datasets)
# ══════════════════════════════════════════════════════════════════════════════

interaction_markers <- c(TARGET_GENE, "MUC5B", "MMP7", "SPP1", "HHIP",
                          "SCGB1A1")

collect_bulk_correlations <- function() {
  results <- list()
  
  # Check each processed dataset
  for (geo_id in c("GSE57148", "GSE76925", "GSE47460", "GSE38974")) {
    # Try bulk VST
    vst_path <- file.path(proc_bulk, geo_id, "vst_matrix.csv")
    # Try array expr
    expr_path <- file.path(proc_array, geo_id, "expr_quantile_gene.csv")
    
    path <- if (file.exists(vst_path)) vst_path else expr_path
    if (!file.exists(path)) next
    
    expr <- read.csv(path, row.names = 1, check.names = FALSE)
    genes_avail <- intersect(interaction_markers, rownames(expr))
    
    if (TARGET_GENE %in% genes_avail && length(genes_avail) >= 2) {
      sub_expr <- t(expr[genes_avail, ])
      
      # Pairwise Spearman with TMEM164
      for (g in setdiff(genes_avail, TARGET_GENE)) {
        ct <- cor.test(sub_expr[, TARGET_GENE], sub_expr[, g],
                       method = "spearman")
        results[[paste0(geo_id, "_", g)]] <- data.frame(
          dataset = geo_id,
          gene    = g,
          rho     = ct$estimate,
          pvalue  = ct$p.value,
          n       = nrow(sub_expr),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results) > 0) rbindlist(results) else NULL
}

cor_df <- collect_bulk_correlations()

if (!is.null(cor_df) && nrow(cor_df) > 0) {
  write.csv(cor_df, file.path(table_dir, "copd_tmem164_marker_correlations.csv"),
            row.names = FALSE)
  
  # Heatmap of correlations
  cor_wide <- dcast(cor_df, gene ~ dataset, value.var = "rho")
  cor_mat <- as.matrix(cor_wide[, -1, drop = FALSE])
  rownames(cor_mat) <- cor_wide$gene
  
  pdf(file.path(fig_dir, "copd_tmem164_marker_correlation_heatmap.pdf"),
      width = 5, height = 3.5)
  pheatmap(cor_mat,
           color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = TRUE, number_format = "%.2f",
           fontsize = 7, main = paste0(TARGET_GENE, " — Spearman ρ with COPD markers"))
  dev.off()
  
  message("✔ Bulk correlation analysis saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. scRNA-seq co-expression (cell-level)
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "COPD_sc_annotated.rds")
if (file.exists(seu_path)) {
  seu <- readRDS(seu_path)
  
  genes_sc <- intersect(interaction_markers, rownames(seu))
  
  if (TARGET_GENE %in% genes_sc && length(genes_sc) >= 2) {
    # Extract expression
    expr_data <- FetchData(seu, vars = c(genes_sc, "celltype_predicted",
                                          "condition"))
    
    # Scatter plots: TMEM164 vs each marker, coloured by condition
    scatter_plots <- list()
    for (g in setdiff(genes_sc, TARGET_GENE)) {
      p <- ggplot(expr_data, aes_string(x = TARGET_GENE, y = g,
                                         colour = "condition")) +
        geom_point(size = 0.1, alpha = 0.3) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
        scale_color_manual(values = col_condition) +
        labs(title = paste0(TARGET_GENE, " vs ", g)) +
        theme_nature(base_size = 6) +
        theme(legend.position = "none")
      scatter_plots[[g]] <- p
    }
    
    if (length(scatter_plots) > 0) {
      p_combined <- wrap_plots(scatter_plots, ncol = 3) +
        plot_annotation(title = paste0(TARGET_GENE,
                                       " co-expression with COPD markers (scRNA-seq)"))
      save_nature_fig(p_combined, "copd_tmem164_coexpression_scatter",
                      width = 7.2, height = 5)
    }
    
    # Per-cell-type correlation heatmap
    ct_cor <- list()
    for (ct in unique(expr_data$celltype_predicted)) {
      sub <- expr_data[expr_data$celltype_predicted == ct, genes_sc]
      if (nrow(sub) > 50) {
        cor_mat <- cor(sub, method = "spearman", use = "complete.obs")
        ct_cor[[ct]] <- cor_mat[TARGET_GENE, ]
      }
    }
    
    if (length(ct_cor) > 0) {
      ct_cor_mat <- do.call(rbind, ct_cor)
      
      pdf(file.path(fig_dir, "copd_tmem164_celltype_correlation.pdf"),
          width = 5, height = 4)
      pheatmap(ct_cor_mat,
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               breaks = seq(-1, 1, length.out = 101),
               display_numbers = TRUE, number_format = "%.2f",
               fontsize = 6,
               main = paste0(TARGET_GENE, " correlation per cell type"))
      dev.off()
    }
  }
  
  rm(seu); gc()
  message("✔ scRNA-seq co-expression analysis saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. WGCNA: identify TMEM164-containing module (bulk data)
# ══════════════════════════════════════════════════════════════════════════════

# Use the largest COPD bulk dataset for WGCNA
wgcna_input <- NULL
for (path in c(file.path(proc_bulk, "GSE57148", "vst_matrix.csv"),
               file.path(proc_array, "GSE76925", "expr_quantile_gene.csv"))) {
  if (file.exists(path)) {
    wgcna_input <- read.csv(path, row.names = 1, check.names = FALSE)
    break
  }
}

if (!is.null(wgcna_input) && requireNamespace("WGCNA", quietly = TRUE)) {
  message("\n── WGCNA module analysis ──")
  
  library(WGCNA)
  allowWGCNAThreads()
  
  # Use top 5000 most variable genes for efficiency
  gene_var <- apply(wgcna_input, 1, var)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]
  if (TARGET_GENE %in% rownames(wgcna_input) &&
      !(TARGET_GENE %in% top_genes)) {
    top_genes <- c(top_genes, TARGET_GENE)
  }
  
  datExpr <- t(wgcna_input[top_genes, ])
  
  # Pick soft threshold
  powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
  power <- sft$powerEstimate
  if (is.na(power)) power <- 6
  
  # Build network
  net <- blockwiseModules(datExpr, power = power,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          verbose = 0)
  
  # Find TMEM164's module
  module_labels <- net$colors
  names(module_labels) <- colnames(datExpr)
  
  if (TARGET_GENE %in% names(module_labels)) {
    tmem_module <- module_labels[TARGET_GENE]
    module_genes <- names(module_labels[module_labels == tmem_module])
    
    message("  TMEM164 module: ", tmem_module,
            " (", length(module_genes), " genes)")
    
    # Save module gene list
    write.csv(data.frame(gene = module_genes, module = tmem_module),
              file.path(table_dir, "copd_tmem164_wgcna_module.csv"),
              row.names = FALSE)
  }
  
  message("  ✔ WGCNA analysis saved")
}

message("\n✔ COPD marker interaction analysis complete")
save_session_info(report_dir)
