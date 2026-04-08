###############################################################################
# celltype_annotation.R — Cell type annotation using canonical markers
#
# Strategy:
#   1. AddModuleScore for each cell lineage (markers from 06/07 files)
#   2. Assign cell type by highest module score
#   3. Validate with DotPlot + FeaturePlot
#   4. Manual review checkpoints
#
# Operates on integrated Seurat objects from preprocess_scrnaseq.R:
#   - data/processed/scrnaseq/COPD_sc_integrated.rds
#   - data/processed/scrnaseq/IPF_sc_integrated.rds
#
# Run: Rscript scripts/02_preprocessing/celltype_annotation.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ══════════════════════════════════════════════════════════════════════════════
# Annotation function
# ══════════════════════════════════════════════════════════════════════════════

annotate_celltypes <- function(seu, obj_name = "dataset") {
  message("\n── Cell type annotation: ", obj_name, " ──")
  message("  Cells: ", ncol(seu), "  Clusters: ",
          length(unique(Idents(seu))))
  
  # --- 1. Module scores for each cell lineage ---
  for (ct in names(markers_all)) {
    # Only use genes present in the dataset
    genes_present <- intersect(markers_all[[ct]], rownames(seu))
    if (length(genes_present) >= 2) {
      seu <- AddModuleScore(seu, features = list(genes_present),
                            name = paste0("score_", ct),
                            seed = 2026)
    } else {
      message("  ⚠ ", ct, ": only ", length(genes_present),
              " markers found — skipping module score")
    }
  }
  
  # --- 2. Assign cell type by highest module score ---
  score_cols <- grep("^score_", colnames(seu@meta.data), value = TRUE)
  if (length(score_cols) > 0) {
    score_mat <- seu@meta.data[, score_cols]
    # Clean column names (remove trailing "1" from AddModuleScore)
    ct_names <- gsub("^score_(.+)1$", "\\1", colnames(score_mat))
    colnames(score_mat) <- ct_names
    
    seu$celltype_predicted <- ct_names[apply(score_mat, 1, which.max)]
    seu$celltype_max_score <- apply(score_mat, 1, max)
    
    # Flag low-confidence annotations
    seu$celltype_confident <- seu$celltype_max_score > 0
    
    tab <- table(seu$celltype_predicted)
    message("  Cell type distribution:")
    for (ct in sort(names(tab))) {
      message("    ", ct, ": ", tab[ct], " (",
              round(100 * tab[ct] / ncol(seu), 1), "%)")
    }
  }
  
  # --- 3. Validation plots ---
  
  # 3a. DotPlot: canonical markers × clusters
  top_markers <- unique(unlist(markers_all))
  top_markers <- intersect(top_markers, rownames(seu))
  
  p_dot <- DotPlot(seu, features = top_markers, group.by = "celltype_predicted",
                   dot.scale = 4) +
    RotatedAxis() +
    labs(title = paste0(obj_name, " — Marker expression by predicted cell type")) +
    theme_nature(base_size = 5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 4))
  
  # 3b. UMAP coloured by cell type
  p_umap <- DimPlot(seu, group.by = "celltype_predicted",
                    cols = col_celltype, pt.size = 0.1, label = TRUE,
                    label.size = 2) +
    labs(title = paste0(obj_name, " — Cell types")) +
    theme_nature() +
    NoLegend()
  
  # 3c. FeaturePlot for TMEM164
  p_tmem164 <- FeaturePlot(seu, features = TARGET_GENE, pt.size = 0.1,
                            order = TRUE) +
    scale_color_viridis_c() +
    labs(title = paste0(TARGET_GENE, " expression")) +
    theme_nature()
  
  # Combined QC panel
  p_combined <- (p_umap | p_tmem164) / p_dot +
    plot_layout(heights = c(1, 2))
  
  save_nature_fig(p_combined,
                  paste0(obj_name, "_celltype_annotation"),
                  width = 7.2, height = 9, dir = report_dir)
  
  message("  ✔ Annotation plots saved")
  return(seu)
}

# ══════════════════════════════════════════════════════════════════════════════
# Annotate COPD integrated object
# ══════════════════════════════════════════════════════════════════════════════

copd_path <- file.path(proc_sc, "COPD_sc_integrated.rds")
if (file.exists(copd_path)) {
  copd_seu <- readRDS(copd_path)
  copd_seu <- annotate_celltypes(copd_seu, "COPD")
  saveRDS(copd_seu, file.path(proc_sc, "COPD_sc_annotated.rds"))
  message("✔ COPD annotated object saved")
} else {
  message("⚠ COPD integrated object not found — run preprocess_scrnaseq.R first")
}

# ══════════════════════════════════════════════════════════════════════════════
# Annotate IPF integrated object
# ══════════════════════════════════════════════════════════════════════════════

ipf_path <- file.path(proc_sc, "IPF_sc_integrated.rds")
if (file.exists(ipf_path)) {
  ipf_seu <- readRDS(ipf_path)
  ipf_seu <- annotate_celltypes(ipf_seu, "IPF")
  saveRDS(ipf_seu, file.path(proc_sc, "IPF_sc_annotated.rds"))
  message("✔ IPF annotated object saved")
} else {
  message("⚠ IPF integrated object not found — run preprocess_scrnaseq.R first")
}

# ══════════════════════════════════════════════════════════════════════════════
# Annotate individual datasets (for per-dataset analyses)
# ══════════════════════════════════════════════════════════════════════════════

individual_rds <- list.files(proc_sc, pattern = "_processed\\.rds$",
                             recursive = TRUE, full.names = TRUE)
for (rds_path in individual_rds) {
  geo_id <- basename(dirname(rds_path))
  message("\n── Individual annotation: ", geo_id, " ──")
  seu <- readRDS(rds_path)
  seu <- tryCatch(
    annotate_celltypes(seu, geo_id),
    error = function(e) {
      message("  ⚠ Annotation failed: ", e$message)
      seu
    }
  )
  saveRDS(seu, sub("_processed\\.rds$", "_annotated.rds", rds_path))
}

# ══════════════════════════════════════════════════════════════════════════════
message("\n✔ Cell type annotation complete")
save_session_info(report_dir)
