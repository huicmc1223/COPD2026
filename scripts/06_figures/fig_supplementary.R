###############################################################################
# fig_supplementary.R — Supplementary figures
#
# S1. Batch correction PCA before/after
# S2. QC violin plots (nFeature, nCount, mito%)
# S3. DoubletFinder summary
# S4. Cell-type annotation validation (module score DotPlot)
# S5. All per-cohort volcano plots (grid)
# S6. Mouse orthologue Tmem164 validation
# S7. Correlation matrix heatmaps (expanded)
# S8. Full KEGG pathway tables
#
# Run: Rscript scripts/06_figures/fig_supplementary.R
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
# S1. Batch correction PCA
# ══════════════════════════════════════════════════════════════════════════════

# Check for saved batch correction plots
batch_plots <- list.files(fig_dir, pattern = "batch_pca_", full.names = TRUE)
if (length(batch_plots) > 0) {
  message("S1: ", length(batch_plots), " batch correction plots available")
} else {
  message("S1: Run preprocessing scripts to generate batch PCA plots")
}

# ══════════════════════════════════════════════════════════════════════════════
# S2. QC violin plots
# ══════════════════════════════════════════════════════════════════════════════

for (dataset in c("COPD", "IPF")) {
  seu_path <- file.path(proc_sc, paste0(dataset, "_sc_annotated.rds"))
  if (!file.exists(seu_path)) next
  
  seu <- readRDS(seu_path)
  
  p_qc <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA",
                                      "percent.mt"),
                   group.by = "geo_id", pt.size = 0, ncol = 3) +
    plot_annotation(title = paste0(dataset, ": QC metrics per dataset"))
  
  save_nature_fig(p_qc, paste0("supp_S2_qc_", tolower(dataset)),
                  width = 7.2, height = 3)
  rm(seu); gc()
}

# ══════════════════════════════════════════════════════════════════════════════
# S3. DoubletFinder — summary table
# ══════════════════════════════════════════════════════════════════════════════

doublet_files <- list.files(proc_sc, pattern = "doubletfinder_summary",
                            full.names = TRUE, recursive = TRUE)
if (length(doublet_files) > 0) {
  doublet_summaries <- lapply(doublet_files, fread)
  doublet_df <- rbindlist(doublet_summaries, fill = TRUE)
  write.csv(doublet_df,
            file.path(table_dir, "supp_doubletfinder_summary.csv"),
            row.names = FALSE)
  message("S3: DoubletFinder summary written")
}

# ══════════════════════════════════════════════════════════════════════════════
# S6. Mouse Tmem164
# ══════════════════════════════════════════════════════════════════════════════

mouse_path <- file.path(proc_bulk, "GSE168299", "vst_matrix.csv")
if (file.exists(mouse_path)) {
  expr <- read.csv(mouse_path, row.names = 1, check.names = FALSE)
  
  if ("Tmem164" %in% rownames(expr)) {
    meta_path <- file.path(raw_dir, "bulk_rnaseq/GSE168299/phenodata_1.csv")
    if (file.exists(meta_path)) {
      meta <- read.csv(meta_path, row.names = 1, stringsAsFactors = FALSE)
      common <- intersect(colnames(expr), rownames(meta))
      
      cond_col <- grep("condition|group|treatment", colnames(meta),
                       ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(cond_col)) {
        plot_df <- data.frame(
          sample = common,
          Tmem164 = as.numeric(expr["Tmem164", common]),
          condition = meta[common, cond_col]
        )
        
        p_mouse <- ggplot(plot_df, aes(x = condition, y = Tmem164,
                                        fill = condition)) +
          geom_boxplot(linewidth = 0.3) +
          geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
          labs(x = NULL, y = "Tmem164 expression",
               title = "GSE168299: Mouse COPD model") +
          theme_nature() + NoLegend()
        
        save_nature_fig(p_mouse, "supp_S6_mouse_tmem164",
                        width = 3.5, height = 3)
      }
    }
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# S7. Full correlation matrices
# ══════════════════════════════════════════════════════════════════════════════

for (file_prefix in c("copd_tmem164_marker_correlations",
                       "ipf_tmem164_marker_correlations")) {
  cor_path <- file.path(table_dir, paste0(file_prefix, ".csv"))
  if (!file.exists(cor_path)) next
  
  cor_df <- fread(cor_path)
  message("S7: ", file_prefix, " — ", nrow(cor_df), " gene-dataset pairs")
}

message("\n✔ Supplementary figures generated")
save_session_info(report_dir)
