###############################################################################
# epithelial_fate.R — Epithelial cell fate & differentiation analysis
#
# Step 23: Signalling axis profiling in COPD/IPF epithelium
#
# Axes: TGFb, Wnt, Notch, Hedgehog, Hippo/YAP, NFkB, mTOR, Ferroptosis,
#        Glycolysis, LipidMetab, Autophagy, ER_Stress
#
# Run: Rscript scripts/05_integrative/epithelial_fate.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Load integrated objects
# ══════════════════════════════════════════════════════════════════════════════

copd_path <- file.path(proc_sc, "COPD_sc_annotated.rds")
ipf_path  <- file.path(proc_sc, "IPF_sc_annotated.rds")

objects <- list()
if (file.exists(copd_path)) objects[["COPD"]] <- readRDS(copd_path)
if (file.exists(ipf_path))  objects[["IPF"]]  <- readRDS(ipf_path)

if (length(objects) == 0) {
  message("⚠ No annotated scRNA-seq objects found")
  quit(save = "no")
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Score epithelial cells for each signalling axis
# ══════════════════════════════════════════════════════════════════════════════

epi_types <- c("AT1", "AT2", "Basal", "Club", "Goblet", "Ciliated")

for (disease in names(objects)) {
  seu <- objects[[disease]]
  epi <- subset(seu, celltype_predicted %in% epi_types)
  
  if (ncol(epi) < 50) next
  message("\n── ", disease, ": ", ncol(epi), " epithelial cells ──")
  
  # Add module scores for each epithelial fate axis
  for (axis_name in names(epithelial_fate)) {
    genes_present <- intersect(epithelial_fate[[axis_name]], rownames(epi))
    if (length(genes_present) >= 3) {
      epi <- AddModuleScore(epi, features = list(genes_present),
                            name = paste0(axis_name, "_"), seed = 2026)
    }
  }
  
  # Collect scores
  score_cols <- grep("_1$", colnames(epi@meta.data), value = TRUE)
  score_cols <- score_cols[grepl(paste(names(epithelial_fate), collapse = "|"),
                                  score_cols)]
  
  if (length(score_cols) == 0) next
  
  # Get TMEM164 expression
  if (TARGET_GENE %in% rownames(epi)) {
    tmem_expr <- FetchData(epi, vars = TARGET_GENE)[, 1]
    median_val <- median(tmem_expr)
    epi$TMEM164_group <- ifelse(tmem_expr > median_val, "High", "Low")
    
    # 2. Heatmap: average score per cell type × axis × TMEM164 group
    agg_list <- list()
    for (sc in score_cols) {
      axis_name <- gsub("_1$", "", sc)
      agg <- aggregate(epi@meta.data[[sc]],
                       by = list(celltype = epi$celltype_predicted,
                                 group = epi$TMEM164_group),
                       FUN = mean)
      colnames(agg)[3] <- "score"
      agg$axis <- axis_name
      agg$disease <- disease
      agg_list[[paste0(disease, "_", sc)]] <- agg
    }
    
    agg_df <- rbindlist(agg_list)
    
    # Heatmap: cell type × axis, split by TMEM164 group
    for (grp in c("High", "Low")) {
      sub <- agg_df[group == grp]
      if (nrow(sub) == 0) next
      
      wide <- dcast(sub, celltype ~ axis, value.var = "score")
      mat <- as.matrix(wide[, -1, drop = FALSE])
      rownames(mat) <- wide$celltype
      
      pdf(file.path(fig_dir,
                    paste0("epithelial_fate_", tolower(disease),
                           "_tmem164_", tolower(grp), ".pdf")),
          width = 6, height = 4)
      pheatmap(mat, scale = "column",
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               fontsize = 6,
               main = paste0(disease, ": ", TARGET_GENE, "-", grp,
                             " epithelial signalling"))
      dev.off()
    }
    
    # 3. Violin plots: top axes with significant difference
    sig_axes <- list()
    for (sc in score_cols) {
      axis_name <- gsub("_1$", "", sc)
      pval <- wilcox.test(epi@meta.data[[sc]] ~ epi$TMEM164_group)$p.value
      if (pval < 0.05) {
        sig_axes[[axis_name]] <- pval
      }
    }
    
    if (length(sig_axes) > 0) {
      sig_df <- data.frame(axis = names(sig_axes),
                           pvalue = unlist(sig_axes))
      sig_df$padj <- p.adjust(sig_df$pvalue, method = "BH")
      write.csv(sig_df,
                file.path(table_dir,
                          paste0("epithelial_fate_", tolower(disease),
                                 "_sig_axes.csv")),
                row.names = FALSE)
      
      message("  Significant axes: ", paste(sig_df$axis, collapse = ", "))
      
      # Violin for top 6 significant axes
      top_axes <- head(sig_df$axis[order(sig_df$pvalue)], 6)
      violin_plots <- list()
      for (ax in top_axes) {
        sc_col <- paste0(ax, "_1")
        p <- VlnPlot(epi, features = sc_col, group.by = "TMEM164_group",
                     cols = c(High = "#E64B35", Low = "#4DBBD5"),
                     pt.size = 0) +
          labs(title = ax, y = "Score") +
          theme_nature(base_size = 6) + NoLegend()
        violin_plots[[ax]] <- p
      }
      
      p_vln <- wrap_plots(violin_plots, ncol = 3)
      save_nature_fig(p_vln,
                      paste0("epithelial_fate_", tolower(disease),
                             "_violin"),
                      width = 7.2, height = 4)
    }
  }
  
  # Save processed epithelial subset
  saveRDS(epi, file.path(proc_sc,
                          paste0(disease, "_epithelial_fate.rds")))
}

message("\n✔ Epithelial fate analysis complete")
save_session_info(report_dir)
