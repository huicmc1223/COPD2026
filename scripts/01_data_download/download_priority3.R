###############################################################################
# download_priority3.R — Download 7 extension/niche datasets
#
# Priority 3 datasets:
#   19. GSE313006 — COPD Spatial transcriptome (Visium, GOLD 3-4)
#   20. GSE168299 — Mouse COPD scRNA-seq (cross-species validation)
#   21. GSE249584 — COPD scRNA-seq Blood (liquid biopsy angle)
#   22. GSE94555  — IPF scRNA-seq (AT2-specific)
#   23. GSE72073  — IPF Array (Chinese cohort)
#   24. GSE11906  — COPD Array (smoking/COPD complex design)
#   25. GSE213001 — IPF Bulk RNA-seq (Australian cohort)
#
# Run: Rscript scripts/01_data_download/download_priority3.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
source(here::here("scripts/utils/geo_helpers.R"))

message("═══════════════════════════════════════════════════════")
message("  Downloading Priority 3 datasets (7 extension)")
message("═══════════════════════════════════════════════════════")

# ── 19. GSE313006 — COPD Spatial (Visium) ────────────────────────────────────
message("\n[19] GSE313006 — Spatial COPD")
download_geo("GSE313006", data_raw, type = "spatial")

# ── 20. GSE168299 — Mouse COPD scRNA-seq ─────────────────────────────────────
message("\n[20] GSE168299 — scRNA-seq Mouse COPD")
download_geo("GSE168299", data_raw, type = "scrnaseq")

# ── 21. GSE249584 — COPD scRNA-seq (Blood) ───────────────────────────────────
message("\n[21] GSE249584 — scRNA-seq COPD Blood")
download_geo("GSE249584", data_raw, type = "scrnaseq")

# ── 22. GSE94555 — IPF scRNA-seq (AT2 cells) ────────────────────────────────
message("\n[22] GSE94555 — scRNA-seq IPF (AT2)")
download_geo("GSE94555", data_raw, type = "scrnaseq")

# ── 23. GSE72073 — IPF Array (Chinese) ──────────────────────────────────────
message("\n[23] GSE72073 — Array IPF (Chinese)")
download_geo("GSE72073", data_raw, type = "microarray")

# ── 24. GSE11906 — COPD Array (smoking design) ──────────────────────────────
message("\n[24] GSE11906 — Array COPD (smoking)")
download_geo("GSE11906", data_raw, type = "microarray")

# ── 25. GSE213001 — IPF Bulk RNA-seq (Australian) ───────────────────────────
message("\n[25] GSE213001 — Bulk RNA-seq IPF (Australian)")
download_geo("GSE213001", data_raw, type = "bulk_rnaseq")

message("\n═══════════════════════════════════════════════════════")
message("  Priority 3 download complete.")
message("  All datasets downloaded. Proceed to preprocessing.")
message("═══════════════════════════════════════════════════════")

save_session_info(report_dir)
