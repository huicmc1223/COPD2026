###############################################################################
# ipf_emt_trajectory.R — EMT trajectory: AT2 → KRT8+ → AT1
#
# Step 17: Pseudotime analysis with TMEM164 expression dynamics
#
# Methods: Monocle3 (primary) or Slingshot (fallback)
#
# Run: Rscript scripts/04_ipf_analysis/ipf_emt_trajectory.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Load IPF annotated object & subset epithelial cells
# ══════════════════════════════════════════════════════════════════════════════

seu_path <- file.path(proc_sc, "IPF_sc_annotated.rds")
if (!file.exists(seu_path)) {
  message("⚠ IPF annotated object not found")
  quit(save = "no")
}

seu <- readRDS(seu_path)

# Subset to epithelial lineages relevant to AT2→AT1 transition
epi_types <- c("AT2", "AT1", "Basal", "Club")
epi <- subset(seu, celltype_predicted %in% epi_types)
message("Epithelial subset: ", ncol(epi), " cells")

# ══════════════════════════════════════════════════════════════════════════════
# EMT module score
# ══════════════════════════════════════════════════════════════════════════════

emt_genes <- intersect(markers_ipf_specific$EMT_markers, rownames(epi))
if (length(emt_genes) >= 3) {
  epi <- AddModuleScore(epi, features = list(emt_genes),
                        name = "EMT_score", seed = 2026)
}

# ══════════════════════════════════════════════════════════════════════════════
# Monocle3 trajectory analysis
# ══════════════════════════════════════════════════════════════════════════════

run_monocle3 <- function(seu_obj) {
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    message("⚠ monocle3 not installed — skipping")
    return(NULL)
  }
  
  library(monocle3)
  library(SeuratWrappers)
  
  message("\n── Monocle3 trajectory analysis ──")
  
  # Convert Seurat → CDS
  cds <- as.cell_data_set(seu_obj)
  
  # Pre-process
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- reduce_dimension(cds, preprocess_method = "PCA")
  
  # Cluster
  cds <- cluster_cells(cds)
  
  # Learn graph
  cds <- learn_graph(cds, use_partition = FALSE)
  
  # Order cells: root = AT2 cells
  at2_cells <- colnames(seu_obj)[seu_obj$celltype_predicted == "AT2"]
  root_cells <- at2_cells[1:min(5, length(at2_cells))]
  cds <- order_cells(cds, root_cells = root_cells)
  
  return(cds)
}

cds <- tryCatch(run_monocle3(epi), error = function(e) {
  message("Monocle3 failed: ", e$message, "\nTrying Slingshot...")
  NULL
})

# ══════════════════════════════════════════════════════════════════════════════
# Slingshot fallback
# ══════════════════════════════════════════════════════════════════════════════

if (is.null(cds) && requireNamespace("slingshot", quietly = TRUE)) {
  library(slingshot)
  message("\n── Slingshot trajectory analysis ──")
  
  # Use UMAP embeddings (or PCA)
  umap_emb <- Embeddings(epi, "umap")
  clusters <- epi$celltype_predicted
  
  sds <- slingshot(umap_emb, clusterLabels = clusters,
                   start.clus = "AT2", end.clus = "AT1")
  
  # Extract pseudotime
  pt <- slingPseudotime(sds)
  epi$pseudotime <- pt[, 1]
  
  message("  Pseudotime range: ", round(min(epi$pseudotime, na.rm = TRUE), 2),
          " — ", round(max(epi$pseudotime, na.rm = TRUE), 2))
}

# ══════════════════════════════════════════════════════════════════════════════
# TMEM164 along trajectory
# ══════════════════════════════════════════════════════════════════════════════

if ("pseudotime" %in% colnames(epi@meta.data) ||
    !is.null(cds)) {
  
  # Extract pseudotime from Monocle3 if available
  if (!is.null(cds)) {
    epi$pseudotime <- monocle3::pseudotime(cds)
  }
  
  # Remove cells with NA pseudotime
  valid <- !is.na(epi$pseudotime)
  epi_pt <- subset(epi, cells = colnames(epi)[valid])
  
  # Get TMEM164 expression
  if (TARGET_GENE %in% rownames(epi_pt)) {
    expr_df <- FetchData(epi_pt, vars = c(TARGET_GENE, "pseudotime",
                                           "celltype_predicted", "condition"))
    if ("EMT_score1" %in% colnames(epi_pt@meta.data)) {
      expr_df$EMT_score <- epi_pt$EMT_score1
    }
    
    # 1. TMEM164 vs pseudotime (scatter + LOESS)
    p_traj <- ggplot(expr_df, aes_string(x = "pseudotime", y = TARGET_GENE)) +
      geom_point(aes(colour = celltype_predicted), size = 0.3, alpha = 0.3) +
      geom_smooth(method = "loess", colour = "#E64B35", linewidth = 1) +
      scale_color_manual(values = col_celltype) +
      labs(x = "Pseudotime (AT2 → AT1)",
           y = paste0(TARGET_GENE, " expression"),
           title = paste0(TARGET_GENE, " along EMT trajectory")) +
      theme_nature()
    
    # 2. EMT score vs pseudotime
    if ("EMT_score" %in% colnames(expr_df)) {
      p_emt <- ggplot(expr_df, aes(x = pseudotime, y = EMT_score)) +
        geom_point(aes(colour = celltype_predicted), size = 0.3, alpha = 0.3) +
        geom_smooth(method = "loess", colour = "#3C5488", linewidth = 1) +
        scale_color_manual(values = col_celltype) +
        labs(x = "Pseudotime", y = "EMT score",
             title = "EMT score along trajectory") +
        theme_nature()
    } else {
      p_emt <- ggplot() + theme_void()
    }
    
    # 3. UMAP coloured by pseudotime
    p_umap_pt <- FeaturePlot(epi_pt, features = "pseudotime", pt.size = 0.3) +
      scale_color_viridis_c(option = "magma") +
      labs(title = "Pseudotime") +
      theme_nature()
    
    # 4. Comparison: IPF vs Control trajectory
    p_traj_cond <- ggplot(expr_df,
                          aes_string(x = "pseudotime", y = TARGET_GENE,
                                     colour = "condition")) +
      geom_smooth(method = "loess", linewidth = 1) +
      scale_color_manual(values = col_condition) +
      labs(x = "Pseudotime", y = paste0(TARGET_GENE, " expression"),
           title = "IPF vs Control trajectory") +
      theme_nature()
    
    # Combined figure
    p_combined <- (p_umap_pt | p_traj) / (p_emt | p_traj_cond)
    save_nature_fig(p_combined, "ipf_emt_trajectory_tmem164",
                    width = 7.2, height = 6)
    
    # Save trajectory data
    write.csv(expr_df, file.path(table_dir, "ipf_emt_trajectory_data.csv"),
              row.names = FALSE)
    
    message("  ✔ EMT trajectory analysis saved")
  }
}

# Save updated object
saveRDS(epi, file.path(proc_sc, "IPF_epithelial_trajectory.rds"))

message("\n✔ IPF EMT trajectory analysis complete")
save_session_info(report_dir)
