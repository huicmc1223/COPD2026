###############################################################################
# install_packages.R — Install all required R packages via renv
# Run once: source(here::here("scripts/00_setup/install_packages.R"))
###############################################################################

# ── Initialize renv (if not yet) ──────────────────────────────────────────────
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::init(bare = TRUE)

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  # Core
  "here", "tidyverse", "data.table", "R.utils",
  # Visualization
  "ggplot2", "patchwork", "cowplot", "ggrepel", "ggsci",
  "pheatmap", "ComplexHeatmap", "circlize", "VennDiagram",
  "ggridges", "ggbeeswarm", "viridis",
  # Statistics
  "metafor", "pROC", "glmnet", "survival", "survminer",
  # Batch correction
  "sva",
  # Microarray
  "lumi",
  # Misc
  "Matrix", "Rcpp", "future", "future.apply", "progressr",
  "writexl", "openxlsx", "clusterProfiler"
)

# ── Bioconductor packages ────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  # GEO data access
  "GEOquery",
  # Bulk RNA-seq
  "DESeq2", "edgeR", "limma",
  # Microarray processing
  "oligo", "affy", "annotate", "AnnotationDbi",
  # Platform annotations
  "hgu133plus2.db", "hugene10stv1.db",
  "illuminaHumanv4.db",
  # Gene sets & pathway analysis
  "org.Hs.eg.db", "org.Mm.eg.db",
  "clusterProfiler", "enrichplot", "DOSE", "fgsea",
  "ReactomePA", "pathview", "msigdbr",
  # Single-cell
  "Seurat", "SeuratObject",
  "harmony",
  "SingleCellExperiment", "scater", "scran",
  "DoubletFinder",
  "AUCell",
  # Trajectory

  "monocle3", "slingshot",
  # Cell communication
  "CellChat", "NicheNet",
  # Network
  "WGCNA", "hdWGCNA",
  # Spatial (Priority 3)
  "Giotto", "SpatialExperiment",
  # Utilities
  "BiocParallel", "GenomicRanges", "SummarizedExperiment",
  "biomaRt", "BPCells"
)

# ── Install ───────────────────────────────────────────────────────────────────
message("Installing CRAN packages...")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

message("Installing Bioconductor packages...")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# ── Snapshot ──────────────────────────────────────────────────────────────────
renv::snapshot()

message("✔ All packages installed. renv snapshot saved.")
