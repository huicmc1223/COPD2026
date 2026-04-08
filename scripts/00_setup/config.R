###############################################################################
# config.R — Global configuration for TMEM164 COPD/IPF analysis project
# Source this file at the top of every script:
#   source(here::here("scripts/00_setup/config.R"))
###############################################################################

# ── Packages ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(here)        # project-relative paths
})

# ── Reproducibility ───────────────────────────────────────────────────────────
set.seed(2026)

# ── Directory paths ───────────────────────────────────────────────────────────
proj_dir      <- here::here()
data_raw      <- here("data", "raw")
data_proc     <- here("data", "processed")
raw_bulk      <- here("data", "raw", "bulk_rnaseq")
raw_array     <- here("data", "raw", "microarray")
raw_sc        <- here("data", "raw", "scrnaseq")
raw_spatial   <- here("data", "raw", "spatial")
raw_meta      <- here("data", "raw", "metadata")
proc_bulk     <- here("data", "processed", "bulk_rnaseq")
proc_array    <- here("data", "processed", "microarray")
proc_sc       <- here("data", "processed", "scrnaseq")
results_dir   <- here("results")
fig_dir       <- here("results", "figures")
table_dir     <- here("results", "tables")
report_dir    <- here("results", "reports")
docs_dir      <- here("docs")
scripts_dir   <- here("scripts")
utils_dir     <- here("scripts", "utils")

# ── Gene of interest ──────────────────────────────────────────────────────────
TARGET_GENE       <- "TMEM164"
TARGET_GENE_MOUSE <- "Tmem164"   # mouse ortholog

# ── Nature-style colour palette ───────────────────────────────────────────────
# Based on NPG (Nature Publishing Group) palette
pal_npg <- c(
  red      = "#E64B35",
  blue     = "#4DBBD5",
  green    = "#00A087",
  navy     = "#3C5488",
  salmon   = "#F39B7F",
  purple   = "#8491B4",
  brown    = "#91D1C2",
  pink     = "#DC0000",
  yellow   = "#7E6148",
  grey     = "#B09C85"
)

# Condition colours (used consistently across all figures)
col_condition <- c(
  Control = "#4DBBD5",
  COPD    = "#E64B35",
  IPF     = "#F39B7F",
  Smoker  = "#8491B4"
)

# GOLD stage colours (sequential: light → dark)
col_gold <- c(
  "GOLD 1" = "#FDD49E",
  "GOLD 2" = "#FDBB84",
  "GOLD 3" = "#FC8D59",
  "GOLD 4" = "#D7301F"
)

# Cell-type colours (23 major types from files 06/07)
col_celltype <- c(
  "AT2"         = "#E64B35", "AT1"         = "#4DBBD5",
  "Basal"       = "#00A087", "Club"        = "#3C5488",
  "Goblet"      = "#F39B7F", "Ciliated"    = "#8491B4",
  "Fibroblast"  = "#91D1C2", "Myofibroblast" = "#DC0000",
  "AM"          = "#7E6148", "IM"          = "#B09C85",
  "Monocyte"    = "#00A087", "DC"          = "#3C5488",
  "CD4 T"       = "#E64B35", "CD8 T"       = "#4DBBD5",
  "B cell"      = "#F39B7F", "NK"          = "#8491B4",
  "Mast"        = "#91D1C2", "Neutrophil"  = "#DC0000",
  "Endothelial" = "#7E6148", "PNEC"        = "#B09C85",
  "Ionocyte"    = "#FDD49E", "Tuft"        = "#FDBB84",
  "Plasma"      = "#FC8D59"
)

# ── Figure dimensions (Nature guidelines) ─────────────────────────────────────
# Single column: 89 mm ≈ 3.5 in; Double column: 183 mm ≈ 7.2 in
fig_single_w  <- 3.5   # inches
fig_double_w  <- 7.2
fig_full_h    <- 9.0   # max height ~ 247 mm ≈ 9.7 in
fig_dpi       <- 300
fig_font_size <- 7     # pt for Nature

# ── Dataset registry ──────────────────────────────────────────────────────────
# Priority 1 — Core datasets
datasets_p1 <- list(
  list(geo_id = "GSE136831", disease = "COPD+IPF", type = "scRNA-seq",
       platform = "10x Chromium", gold = NA,
       n_ctrl = 28, n_case = 50, race = "Caucasian"),
  list(geo_id = "GSE57148",  disease = "COPD", type = "bulk_rnaseq",
       platform = "Illumina HiSeq", gold = "1-2",
       n_ctrl = 91, n_case = 98, race = "Korean"),
  list(geo_id = "GSE76925",  disease = "COPD", type = "microarray",
       platform = "Illumina HT-12 V4", gold = "3-4",
       n_ctrl = 40, n_case = 111, race = "Caucasian"),
  list(geo_id = "GSE310058", disease = "COPD", type = "snRNA-seq",
       platform = "10x Chromium", gold = "1-4",
       n_ctrl = 28, n_case = 113, race = "Mixed"),
  list(geo_id = "GSE47460",  disease = "COPD+IPF", type = "microarray",
       platform = "Agilent 4x44K/8x60K", gold = NA,
       n_ctrl = 17, n_case = 264, race = "Caucasian"),
  list(geo_id = "GSE135893", disease = "IPF", type = "scRNA-seq",
       platform = "10x Chromium", gold = NA,
       n_ctrl = 9, n_case = 21, race = "European"),
  list(geo_id = "GSE134692", disease = "IPF", type = "bulk_rnaseq",
       platform = "Illumina HiSeq", gold = NA,
       n_ctrl = 26, n_case = 46, race = "Caucasian"),
  list(geo_id = "GSE124685", disease = "IPF", type = "bulk_rnaseq",
       platform = "Illumina HiSeq", gold = NA,
       n_ctrl = 35, n_case = 49, race = "Caucasian"),
  list(geo_id = "GSE32537",  disease = "IPF", type = "microarray",
       platform = "Affymetrix HuGene-1_0-st", gold = NA,
       n_ctrl = 50, n_case = 169, race = "Caucasian")
)

# Priority 2 — Validation datasets
datasets_p2 <- list(
  list(geo_id = "GSE38974",  disease = "COPD", type = "microarray",
       platform = "Agilent 4x44K", gold = "1,2,4", race = "Caucasian"),
  list(geo_id = "GSE227691", disease = "COPD", type = "scRNA-seq",
       platform = "10x Chromium", gold = "1,2", race = "Chinese"),
  list(geo_id = "GSE173896", disease = "COPD", type = "scRNA-seq",
       platform = "10x Chromium", gold = "3-4", race = "Japanese"),
  list(geo_id = "GSE279570", disease = "COPD", type = "scRNA-seq",
       platform = "10x Chromium", gold = "2", race = "Chinese"),
  list(geo_id = "GSE171541", disease = "COPD", type = "scRNA-seq",
       platform = "10x Chromium", gold = "1-2", race = "Caucasian"),
  list(geo_id = "GSE128033", disease = "IPF", type = "scRNA-seq",
       platform = "10x Chromium", gold = NA, race = "Caucasian"),
  list(geo_id = "GSE122960", disease = "IPF", type = "scRNA-seq",
       platform = "10x Chromium", gold = NA, race = "Caucasian"),
  list(geo_id = "GSE24206",  disease = "IPF", type = "microarray",
       platform = "Affymetrix HG-U133_Plus_2", gold = NA, race = "Caucasian"),
  list(geo_id = "GSE110147", disease = "IPF", type = "microarray",
       platform = "Affymetrix HuGene-1_0-st", gold = NA, race = "European")
)

# Priority 3 — Extension datasets
datasets_p3 <- list(
  list(geo_id = "GSE313006", disease = "COPD", type = "spatial",
       platform = "Visium", gold = "3-4", race = "Mixed"),
  list(geo_id = "GSE168299", disease = "COPD", type = "scRNA-seq",
       platform = "10x Chromium", gold = NA, race = "Mouse",
       species = "Mus musculus"),
  list(geo_id = "GSE249584", disease = "COPD", type = "scRNA-seq",
       platform = "10x Chromium", gold = "2", race = "Korean",
       tissue = "Blood"),
  list(geo_id = "GSE94555",  disease = "IPF", type = "scRNA-seq",
       platform = "10x Chromium", gold = NA, race = "Caucasian"),
  list(geo_id = "GSE72073",  disease = "IPF", type = "microarray",
       platform = "Affymetrix HTA-2_0", gold = NA, race = "Chinese"),
  list(geo_id = "GSE11906",  disease = "COPD", type = "microarray",
       platform = "Affymetrix HG-U133_Plus_2", gold = "1-2", race = "Mixed"),
  list(geo_id = "GSE213001", disease = "IPF", type = "bulk_rnaseq",
       platform = "Illumina", gold = NA, race = "Australian")
)

# ── QC thresholds (scRNA-seq) ─────────────────────────────────────────────────
sc_min_features <- 200
sc_max_features <- 8000
sc_max_mito_pct <- 20
sc_min_cells    <- 3

# ── Analysis parameters ───────────────────────────────────────────────────────
deseq2_padj_cutoff <- 0.05
deseq2_lfc_cutoff  <- 0.5     # log2 fold change
limma_padj_cutoff  <- 0.05
limma_lfc_cutoff   <- 0.5
gsea_pval_cutoff   <- 0.05

# ── Helper: load all utility scripts ──────────────────────────────────────────
load_utils <- function() {
  util_files <- list.files(utils_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in util_files) source(f)
}

# ── Session info helper ───────────────────────────────────────────────────────
save_session_info <- function(output_dir) {
  si <- sessionInfo()
  writeLines(capture.output(print(si)),
             file.path(output_dir, paste0("session_info_",
                                          format(Sys.time(), "%Y%m%d_%H%M%S"),
                                          ".txt")))
}

message("✔ config.R loaded — project root: ", proj_dir)
