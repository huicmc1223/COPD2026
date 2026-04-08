###############################################################################
# preprocess_scrnaseq.R — scRNA-seq / snRNA-seq preprocessing + integration
#
# Per-dataset pipeline (Seurat v5):
#   CreateSeuratObject → QC → DoubletFinder → SCTransform →
#   PCA → UMAP → Find Clusters → Save .rds
#
# Multi-dataset integration (Harmony):
#   COPD integration: GSE136831(COPD) + GSE310058 + GSE171541 + GSE227691 + GSE173896
#   IPF integration:  GSE136831(IPF) + GSE135893 + GSE128033 + GSE122960
#
# Run: Rscript scripts/02_preprocessing/preprocess_scrnaseq.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

# ══════════════════════════════════════════════════════════════════════════════
# Helper: Standard per-dataset QC + processing
# ══════════════════════════════════════════════════════════════════════════════

preprocess_one_sc <- function(data_path, geo_id, project_name = NULL,
                              min_features = sc_min_features,
                              max_features = sc_max_features,
                              max_mito = sc_max_mito_pct) {
  message("\n── Processing ", geo_id, " ──")
  
  if (is.null(project_name)) project_name <- geo_id
  
  # Read data (10x format or h5)
  h5_files <- list.files(data_path, pattern = "\\.h5$", full.names = TRUE)
  
  if (length(h5_files) > 0) {
    counts <- Read10X_h5(h5_files[1])
    if (is.list(counts)) counts <- counts[["Gene Expression"]]
  } else {
    # Try standard 10x directory structure
    counts <- Read10X(data_path)
  }
  
  seu <- CreateSeuratObject(counts = counts, project = project_name,
                            min.cells = sc_min_cells,
                            min.features = min_features)
  
  # Add percent mitochondrial
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  message("  Raw: ", ncol(seu), " cells × ", nrow(seu), " genes")
  
  # ── QC filtering ──
  seu <- subset(seu,
                subset = nFeature_RNA > min_features &
                  nFeature_RNA < max_features &
                  percent.mt < max_mito)
  
  message("  After QC: ", ncol(seu), " cells")
  
  # ── QC violin plot ──
  p_qc <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0) &
    theme_nature()
  
  # ── Doublet detection (DoubletFinder) ──
  tryCatch({
    if (requireNamespace("DoubletFinder", quietly = TRUE) && ncol(seu) > 100) {
      seu <- NormalizeData(seu, verbose = FALSE)
      seu <- FindVariableFeatures(seu, verbose = FALSE)
      seu <- ScaleData(seu, verbose = FALSE)
      seu <- RunPCA(seu, verbose = FALSE)
      
      # Estimate pK
      sweep.res <- DoubletFinder::paramSweep(seu, PCs = 1:20, sct = FALSE)
      sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
      bcmvn <- DoubletFinder::find.pK(sweep.stats)
      pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      
      # Expected doublet rate (~0.8% per 1000 cells)
      nExp <- round(0.008 * ncol(seu)^2 / 10000)
      
      seu <- DoubletFinder::doubletFinder(seu, PCs = 1:20, pN = 0.25,
                                           pK = pK, nExp = nExp, sct = FALSE)
      
      # Remove doublets
      df_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
      if (length(df_col) > 0) {
        n_doublets <- sum(seu@meta.data[[df_col[1]]] == "Doublet")
        seu <- subset(seu, cells = colnames(seu)[
          seu@meta.data[[df_col[1]]] == "Singlet"])
        message("  Removed ", n_doublets, " doublets → ", ncol(seu), " cells")
      }
    }
  }, error = function(e) {
    message("  ⚠ DoubletFinder skipped: ", e$message)
  })
  
  # ── SCTransform + dimensional reduction ──
  seu <- SCTransform(seu, verbose = FALSE, seed.use = 2026)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
  
  # Add dataset label
  seu$dataset <- geo_id
  
  message("  Clusters: ", length(unique(Idents(seu))))
  
  # ── Save ──
  out_dir <- file.path(proc_sc, geo_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(seu, file.path(out_dir, paste0(geo_id, "_processed.rds")))
  save_nature_fig(p_qc, paste0(geo_id, "_qc_violin"), width = 7, height = 3,
                  dir = report_dir)
  
  message("  ✔ Saved: ", geo_id, "_processed.rds")
  return(seu)
}

# ══════════════════════════════════════════════════════════════════════════════
# Process Priority 1 scRNA-seq datasets
# ══════════════════════════════════════════════════════════════════════════════

sc_objects <- list()

# Dataset paths (adapt after actual download)
sc_datasets <- list(
  GSE136831 = file.path(raw_sc, "GSE136831"),
  GSE310058 = file.path(raw_sc, "GSE310058"),
  GSE135893 = file.path(raw_sc, "GSE135893")
)

for (geo_id in names(sc_datasets)) {
  data_path <- sc_datasets[[geo_id]]
  if (dir.exists(data_path) &&
      length(list.files(data_path, recursive = TRUE)) > 0) {
    sc_objects[[geo_id]] <- tryCatch(
      preprocess_one_sc(data_path, geo_id),
      error = function(e) {
        message("⚠ ", geo_id, " processing failed: ", e$message)
        NULL
      }
    )
  } else {
    message("⚠ ", geo_id, " data not found at ", data_path, " — skipping")
  }
}

# Priority 2 scRNA-seq datasets
p2_sc <- c("GSE227691", "GSE173896", "GSE279570", "GSE171541",
           "GSE128033", "GSE122960")
for (geo_id in p2_sc) {
  data_path <- file.path(raw_sc, geo_id)
  if (dir.exists(data_path) &&
      length(list.files(data_path, recursive = TRUE)) > 0) {
    sc_objects[[geo_id]] <- tryCatch(
      preprocess_one_sc(data_path, geo_id),
      error = function(e) {
        message("⚠ ", geo_id, " failed: ", e$message)
        NULL
      }
    )
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Multi-dataset Harmony integration
# ══════════════════════════════════════════════════════════════════════════════

integrate_datasets <- function(seu_list, integration_name, batch_var = "dataset") {
  message("\n── Harmony integration: ", integration_name, " ──")
  
  # Filter NULLs
  seu_list <- seu_list[!sapply(seu_list, is.null)]
  if (length(seu_list) < 2) {
    message("  ⚠ Need >= 2 datasets for integration; skipping")
    return(NULL)
  }
  
  message("  Merging ", length(seu_list), " datasets...")
  
  # Merge Seurat objects
  merged <- merge(seu_list[[1]], y = seu_list[-1],
                  add.cell.ids = names(seu_list))
  
  message("  Merged: ", ncol(merged), " cells × ", nrow(merged), " genes")
  
  # Re-run SCTransform on merged
  merged <- SCTransform(merged, verbose = FALSE, seed.use = 2026)
  merged <- RunPCA(merged, verbose = FALSE)
  
  # Harmony
  merged <- harmony_integrate(merged, batch_var = batch_var, dims = 30)
  
  # Save
  out_path <- file.path(proc_sc, paste0(integration_name, "_integrated.rds"))
  saveRDS(merged, out_path)
  
  # UMAP plots
  p1 <- DimPlot(merged, group.by = batch_var, pt.size = 0.1) +
    labs(title = paste0(integration_name, " — by Dataset")) +
    theme_nature()
  p2 <- DimPlot(merged, group.by = "seurat_clusters", pt.size = 0.1,
                label = TRUE) +
    labs(title = paste0(integration_name, " — Clusters")) +
    theme_nature()
  
  p_combined <- p1 + p2
  save_nature_fig(p_combined, paste0(integration_name, "_integration_umap"),
                  width = 7, height = 3.5, dir = report_dir)
  
  message("  ✔ Integration saved: ", out_path)
  return(merged)
}

# ── COPD integration ──
# GSE136831 needs to be subset to COPD + Control samples only
copd_sc_datasets <- c("GSE136831", "GSE310058", "GSE171541",
                       "GSE227691", "GSE173896")
copd_sc_list <- sc_objects[intersect(names(sc_objects), copd_sc_datasets)]
if (length(copd_sc_list) >= 2) {
  copd_integrated <- integrate_datasets(copd_sc_list, "COPD_sc")
}

# ── IPF integration ──
ipf_sc_datasets <- c("GSE136831", "GSE135893", "GSE128033", "GSE122960")
ipf_sc_list <- sc_objects[intersect(names(sc_objects), ipf_sc_datasets)]
if (length(ipf_sc_list) >= 2) {
  ipf_integrated <- integrate_datasets(ipf_sc_list, "IPF_sc")
}

# ══════════════════════════════════════════════════════════════════════════════
message("\n✔ scRNA-seq preprocessing complete")
save_session_info(report_dir)
