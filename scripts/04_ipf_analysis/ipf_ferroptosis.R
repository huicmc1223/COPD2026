###############################################################################
# ipf_ferroptosis.R — Ferroptosis, autophagy, lipid metabolism, Hippo/YAP
#
# Step 18: TMEM164 association with key IPF pathogenic pathways
#
# Run: Rscript scripts/04_ipf_analysis/ipf_ferroptosis.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# 1. Bulk-level: TMEM164 correlation with pathway genes
# ══════════════════════════════════════════════════════════════════════════════

pathway_gene_lists <- list(
  Ferroptosis  = ferroptosis_genes,
  Autophagy    = autophagy_genes,
  Hippo_YAP    = hippo_yap_genes,
  Lipid_metab  = c("SREBF1", "PPARG", "LPCAT3", "PLA2G4A", "ACSL4"),
  TGFb_EMT     = epithelial_fate$TGFb_axis
)

bulk_pathway_cor <- function() {
  results <- list()
  
  # Use available IPF bulk datasets
  for (geo_id in c("GSE134692", "GSE124685")) {
    vst_path <- file.path(proc_bulk, geo_id, "vst_matrix.csv")
    if (!file.exists(vst_path)) next
    
    expr <- read.csv(vst_path, row.names = 1, check.names = FALSE)
    if (!(TARGET_GENE %in% rownames(expr))) next
    
    tmem_expr <- as.numeric(expr[TARGET_GENE, ])
    
    for (pw in names(pathway_gene_lists)) {
      pw_genes <- intersect(pathway_gene_lists[[pw]], rownames(expr))
      for (g in pw_genes) {
        ct <- cor.test(tmem_expr, as.numeric(expr[g, ]),
                       method = "spearman")
        results[[paste0(geo_id, "_", pw, "_", g)]] <- data.frame(
          dataset  = geo_id, pathway = pw, gene = g,
          rho      = ct$estimate, pvalue = ct$p.value
        )
      }
    }
  }
  
  if (length(results) > 0) rbindlist(results) else NULL
}

bulk_cor <- bulk_pathway_cor()

if (!is.null(bulk_cor) && nrow(bulk_cor) > 0) {
  bulk_cor$padj <- p.adjust(bulk_cor$pvalue, method = "BH")
  write.csv(bulk_cor,
            file.path(table_dir, "ipf_tmem164_pathway_correlations.csv"),
            row.names = FALSE)
  
  # Heatmap: rows = genes, columns = datasets, values = rho
  for (pw in unique(bulk_cor$pathway)) {
    sub <- bulk_cor[pathway == pw]
    if (nrow(sub) > 0) {
      mat <- dcast(sub, gene ~ dataset, value.var = "rho")
      mat_m <- as.matrix(mat[, -1, drop = FALSE])
      rownames(mat_m) <- mat$gene
      
      pdf(file.path(fig_dir, paste0("ipf_tmem164_", tolower(pw), "_cor.pdf")),
          width = 4, height = max(2.5, nrow(mat_m) * 0.25 + 1))
      pheatmap(mat_m,
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               breaks = seq(-1, 1, length.out = 101),
               display_numbers = TRUE, number_format = "%.2f",
               fontsize = 6,
               main = paste0(TARGET_GENE, " vs ", pw, " genes (Spearman ρ)"))
      dev.off()
    }
  }
  message("✔ Bulk pathway correlations saved")
}

# ══════════════════════════════════════════════════════════════════════════════
# 2. scRNA-seq: Pathway module scores + TMEM164 correlation
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "IPF_sc_annotated.rds")
if (file.exists(seu_path)) {
  library(Seurat)
  seu <- readRDS(seu_path)
  
  # Add module scores for each pathway
  for (pw in names(pathway_gene_lists)) {
    genes_present <- intersect(pathway_gene_lists[[pw]], rownames(seu))
    if (length(genes_present) >= 3) {
      seu <- AddModuleScore(seu, features = list(genes_present),
                            name = paste0(pw, "_score"), seed = 2026)
    }
  }
  
  # Extract TMEM164 + pathway scores
  score_cols <- grep("_score1$", colnames(seu@meta.data), value = TRUE)
  if (TARGET_GENE %in% rownames(seu) && length(score_cols) > 0) {
    plot_data <- FetchData(seu, vars = c(TARGET_GENE, score_cols,
                                          "celltype_predicted", "condition"))
    
    # Scatter: TMEM164 vs each pathway score
    scatter_plots <- list()
    for (sc in score_cols) {
      pw_name <- gsub("_score1$", "", sc)
      p <- ggplot(plot_data, aes_string(x = TARGET_GENE, y = sc)) +
        geom_point(aes(colour = condition), size = 0.1, alpha = 0.2) +
        geom_smooth(method = "lm", linewidth = 0.5, colour = "black") +
        scale_color_manual(values = col_condition) +
        labs(x = TARGET_GENE, y = paste0(pw_name, " score"),
             title = pw_name) +
        theme_nature(base_size = 6) +
        theme(legend.position = "none")
      scatter_plots[[pw_name]] <- p
    }
    
    p_combined <- wrap_plots(scatter_plots, ncol = 3)
    save_nature_fig(p_combined, "ipf_tmem164_pathway_scatter",
                    width = 7.2, height = 5)
    
    # Per-cell-type: mean pathway score in TMEM164-high vs TMEM164-low
    median_expr <- median(plot_data[[TARGET_GENE]])
    plot_data$TMEM164_group <- ifelse(plot_data[[TARGET_GENE]] > median_expr,
                                       "High", "Low")
    
    summary_list <- list()
    for (sc in score_cols) {
      pw_name <- gsub("_score1$", "", sc)
      agg <- aggregate(plot_data[[sc]],
                       by = list(celltype = plot_data$celltype_predicted,
                                 TMEM164_group = plot_data$TMEM164_group,
                                 condition = plot_data$condition),
                       FUN = mean)
      colnames(agg)[4] <- "mean_score"
      agg$pathway <- pw_name
      summary_list[[pw_name]] <- agg
    }
    
    summary_df <- rbindlist(summary_list)
    write.csv(summary_df,
              file.path(table_dir, "ipf_tmem164_pathway_celltype_summary.csv"),
              row.names = FALSE)
    
    # Grouped bar plot
    p_bar <- ggplot(summary_df[condition == "IPF"],
                    aes(x = celltype, y = mean_score,
                        fill = TMEM164_group)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ pathway, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = c(High = "#E64B35", Low = "#4DBBD5")) +
      labs(title = paste0(TARGET_GENE, "-high vs -low: pathway activity in IPF"),
           x = "Cell type", y = "Mean pathway score") +
      theme_nature(base_size = 5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_nature_fig(p_bar, "ipf_tmem164_pathway_barplot",
                    width = 7.2, height = 6)
  }
  
  rm(seu); gc()
  message("✔ scRNA-seq pathway analysis saved")
}

message("\n✔ IPF ferroptosis/autophagy/Hippo analysis complete")
save_session_info(report_dir)
