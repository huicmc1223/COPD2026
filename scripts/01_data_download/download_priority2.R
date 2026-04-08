###############################################################################
# download_priority2.R — Download 9 validation datasets
#
# Priority 2 datasets:
#   10. GSE38974  — COPD Array (GOLD 1,2,4)
#   11. GSE227691 — COPD scRNA-seq (Chinese, GOLD 1-2)
#   12. GSE173896 — COPD scRNA-seq (Japanese, GOLD 3-4)
#   13. GSE279570 — COPD scRNA-seq (Chinese, GOLD 2)
#   14. GSE171541 — COPD scRNA-seq (age-stratified, GOLD 1-2)
#   15. GSE128033 — IPF scRNA-seq
#   16. GSE122960 — IPF scRNA-seq
#   17. GSE24206  — IPF Array (Early vs Advanced)
#   18. GSE110147 — IPF Array (UIP/NSIP)
#
# Run: Rscript scripts/01_data_download/download_priority2.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
source(here::here("scripts/utils/geo_helpers.R"))

message("═══════════════════════════════════════════════════════")
message("  Downloading Priority 2 datasets (9 validation)")
message("═══════════════════════════════════════════════════════")

# ── 10. GSE38974 — COPD Microarray (Agilent, GOLD 1,2,4) ────────────────────
message("\n[10] GSE38974 — Array COPD (multi-GOLD)")
download_geo("GSE38974", data_raw, type = "microarray")

# ── 11. GSE227691 — COPD scRNA-seq (Chinese, mild/moderate) ──────────────────
message("\n[11] GSE227691 — scRNA-seq COPD (Chinese)")
download_geo("GSE227691", data_raw, type = "scrnaseq")

# ── 12. GSE173896 — COPD scRNA-seq (Japanese, GOLD 3-4) ─────────────────────
message("\n[12] GSE173896 — scRNA-seq COPD (Japanese)")
download_geo("GSE173896", data_raw, type = "scrnaseq")

# ── 13. GSE279570 — COPD scRNA-seq (Chinese, GOLD 2) ────────────────────────
message("\n[13] GSE279570 — scRNA-seq COPD (Chinese)")
download_geo("GSE279570", data_raw, type = "scrnaseq")

# ── 14. GSE171541 — COPD scRNA-seq (age-stratified) ─────────────────────────
message("\n[14] GSE171541 — scRNA-seq COPD (age groups)")
download_geo("GSE171541", data_raw, type = "scrnaseq")

# ── 15. GSE128033 — IPF scRNA-seq ────────────────────────────────────────────
message("\n[15] GSE128033 — scRNA-seq IPF")
download_geo("GSE128033", data_raw, type = "scrnaseq")

# ── 16. GSE122960 — IPF scRNA-seq ────────────────────────────────────────────
message("\n[16] GSE122960 — scRNA-seq IPF")
download_geo("GSE122960", data_raw, type = "scrnaseq")

# ── 17. GSE24206 — IPF array (Early vs Advanced) ────────────────────────────
message("\n[17] GSE24206 — Array IPF (staged)")
download_geo("GSE24206", data_raw, type = "microarray")

# ── 18. GSE110147 — IPF array (UIP/NSIP subtypes) ───────────────────────────
message("\n[18] GSE110147 — Array IPF (UIP/NSIP)")
download_geo("GSE110147", data_raw, type = "microarray")

message("\n═══════════════════════════════════════════════════════")
message("  Priority 2 download complete.")
message("═══════════════════════════════════════════════════════")

save_session_info(report_dir)
