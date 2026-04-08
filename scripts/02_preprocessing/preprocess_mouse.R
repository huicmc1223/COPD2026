###############################################################################
# preprocess_mouse.R — GSE168299 Mouse COPD scRNA-seq preprocessing
#
# Cross-species validation: check Tmem164 (mouse ortholog) in murine COPD.
# Pipeline: Standard Seurat → cell type annotation with mouse markers
#
# Run: Rscript scripts/02_preprocessing/preprocess_mouse.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
load_utils()

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ══════════════════════════════════════════════════════════════════════════════
# Mouse orthologs of key markers (capitalise first letter)
# ══════════════════════════════════════════════════════════════════════════════

to_mouse <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

markers_mouse <- lapply(markers_all, function(genes) to_mouse(genes))

# ══════════════════════════════════════════════════════════════════════════════

geo_id <- "GSE168299"
data_path <- file.path(raw_sc, geo_id)

if (!dir.exists(data_path) ||
    length(list.files(data_path, recursive = TRUE)) == 0) {
  message("⚠ GSE168299 data not found — run download_priority3.R first")
  quit(save = "no")
}

message("── Processing GSE168299 — Mouse COPD scRNA-seq ──")

# ── Read data ──
h5_files <- list.files(data_path, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) > 0) {
  counts <- Read10X_h5(h5_files[1])
} else {
  counts <- Read10X(data_path)
}

seu <- CreateSeuratObject(counts = counts, project = "GSE168299_mouse",
                          min.cells = 3, min.features = 200)

# Mouse mitochondrial genes start with "mt-"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

message("  Raw: ", ncol(seu), " cells")

# QC filter
seu <- subset(seu,
              subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                percent.mt < 15)
message("  After QC: ", ncol(seu), " cells")

# Standard pipeline
seu <- SCTransform(seu, verbose = FALSE, seed.use = 2026)
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

# ── Assign condition from metadata ──
# GSE168299: Ctr = GSM5136210-GSM5136213, COPD = GSM5136206-GSM5136209
# This needs to be mapped from downloadeed phenodata — placeholder:
seu$condition <- ifelse(grepl("Ctrl|Control|Ctr", seu$orig.ident,
                              ignore.case = TRUE), "Control", "COPD")

# ── Cell type annotation (mouse markers) ──
for (ct in names(markers_mouse)) {
  genes_present <- intersect(markers_mouse[[ct]], rownames(seu))
  if (length(genes_present) >= 2) {
    seu <- AddModuleScore(seu, features = list(genes_present),
                          name = paste0("score_", ct), seed = 2026)
  }
}

score_cols <- grep("^score_", colnames(seu@meta.data), value = TRUE)
if (length(score_cols) > 0) {
  score_mat <- seu@meta.data[, score_cols]
  ct_names <- gsub("^score_(.+)1$", "\\1", colnames(score_mat))
  seu$celltype_predicted <- ct_names[apply(score_mat, 1, which.max)]
}

# ── Tmem164 analysis ──
p1 <- DimPlot(seu, group.by = "celltype_predicted", pt.size = 0.3,
              label = TRUE, label.size = 2) +
  labs(title = "Mouse COPD — Cell types") + theme_nature() + NoLegend()

p2 <- FeaturePlot(seu, features = TARGET_GENE_MOUSE, pt.size = 0.3,
                  order = TRUE) +
  scale_color_viridis_c() +
  labs(title = paste0(TARGET_GENE_MOUSE, " expression")) +
  theme_nature()

p3 <- VlnPlot(seu, features = TARGET_GENE_MOUSE, group.by = "celltype_predicted",
              split.by = "condition", pt.size = 0) +
  scale_fill_manual(values = col_condition) +
  labs(title = paste0(TARGET_GENE_MOUSE, " per cell type")) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_combined <- (p1 | p2) / p3 + plot_layout(heights = c(1, 1))
save_nature_fig(p_combined, "GSE168299_mouse_tmem164",
                width = 7.2, height = 6, dir = report_dir)

# ── Save ──
out_dir <- file.path(proc_sc, geo_id)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(seu, file.path(out_dir, "GSE168299_mouse_annotated.rds"))

message("✔ Mouse COPD preprocessing complete")
save_session_info(report_dir)
