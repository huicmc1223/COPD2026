###############################################################################
# download_priority1.R — Download 9 core datasets
#
# Priority 1 datasets:
#   1. GSE136831 — COPD+IPF scRNA-seq (28C + 18COPD + 32IPF)
#   2. GSE57148  — COPD Bulk RNA-seq (91C + 98COPD, GOLD 1-2)
#   3. GSE76925  — COPD Array/Illumina (40C + 111COPD, GOLD 3-4)
#   4. GSE310058 — COPD snRNA-seq (141 samples, GOLD 1-4)
#   5. GSE47460  — COPD+IPF Array/Agilent (17C + 75COPD + 189IPF)
#   6. GSE135893 — IPF scRNA-seq (9C + 21IPF)
#   7. GSE134692 — IPF Bulk RNA-seq (26C + 46IPF)
#   8. GSE124685 — IPF Bulk RNA-seq (35C + 49IPF)
#   9. GSE32537  — IPF Array/Affy (50C + 169IPF)
#
# Run: Rscript scripts/01_data_download/download_priority1.R
###############################################################################

source(here::here("scripts/00_setup/config.R"))
source(here::here("scripts/utils/geo_helpers.R"))

message("═══════════════════════════════════════════════════════")
message("  Downloading Priority 1 datasets (9 core)")
message("═══════════════════════════════════════════════════════")

# ── 1. GSE136831 — scRNA-seq (COPD + IPF + Control) ──────────────────────────
# This is a large 10x dataset; supplementary files contain raw counts
message("\n[1/9] GSE136831 — scRNA-seq COPD+IPF")
download_geo("GSE136831", data_raw, type = "scrnaseq")

# ── 2. GSE57148 — Bulk RNA-seq COPD (Korean, GOLD 1-2) ───────────────────────
message("\n[2/9] GSE57148 — Bulk RNA-seq COPD")
download_geo("GSE57148", data_raw, type = "bulk_rnaseq")

# ── 3. GSE76925 — Microarray COPD (Illumina HT-12, GOLD 3-4) ────────────────
message("\n[3/9] GSE76925 — Array COPD")
download_geo("GSE76925", data_raw, type = "microarray")

# ── 4. GSE310058 — snRNA-seq COPD (GOLD 1-4, multi-ethnic) ───────────────────
# Part of GSE313090 super-series
message("\n[4/9] GSE310058 — snRNA-seq COPD")
download_geo("GSE310058", data_raw, type = "scrnaseq")

# ── 5. GSE47460 — Microarray COPD+IPF (Agilent, largest IPF array) ───────────
# NOTE: This dataset has two platforms (GPL6480, GPL14550). Both will download.
message("\n[5/9] GSE47460 — Array COPD+IPF")
download_geo("GSE47460", data_raw, type = "microarray")

# ── 6. GSE135893 — scRNA-seq IPF ──────────────────────────────────────────────
message("\n[6/9] GSE135893 — scRNA-seq IPF")
download_geo("GSE135893", data_raw, type = "scrnaseq")

# ── 7. GSE134692 — Bulk RNA-seq IPF ──────────────────────────────────────────
message("\n[7/9] GSE134692 — Bulk RNA-seq IPF")
download_geo("GSE134692", data_raw, type = "bulk_rnaseq")

# ── 8. GSE124685 — Bulk RNA-seq IPF (validation) ─────────────────────────────
message("\n[8/9] GSE124685 — Bulk RNA-seq IPF")
download_geo("GSE124685", data_raw, type = "bulk_rnaseq")

# ── 9. GSE32537 — Microarray IPF (Affy HuGene-1_0-st) ───────────────────────
message("\n[9/9] GSE32537 — Array IPF")
download_geo("GSE32537", data_raw, type = "microarray")

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n═══════════════════════════════════════════════════════")
message("  Priority 1 download complete.")
message("  Check data/raw/ for downloaded files.")
message("  Next: run download_priority2.R or proceed to preprocessing.")
message("═══════════════════════════════════════════════════════")

save_session_info(report_dir)
