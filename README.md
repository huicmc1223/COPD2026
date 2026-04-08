# COPD2026 — TMEM164 in COPD & IPF

> **TMEM164 as an early diagnostic biomarker for COPD and its role in IPF fibrosis: a multi-cohort transcriptomic study**

[![R](https://img.shields.io/badge/R-≥4.3-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This project systematically investigates **TMEM164** across 25 GEO datasets (18 COPD + 16 IPF, with overlap) spanning bulk RNA-seq, microarray, scRNA-seq, snRNA-seq, and spatial transcriptomics platforms. The analysis follows the 2026 GOLD guidelines (ABE classification) and produces Nature publication-grade figures.

### Key questions
1. Is TMEM164 consistently differentially expressed in COPD vs Control across cohorts?
2. Does TMEM164 expression correlate with COPD severity (GOLD 1→4)?
3. Which cell types drive TMEM164 dysregulation?
4. What is TMEM164's role in the AT2→KRT8+→AT1 epithelial-mesenchymal transition in IPF?
5. Is TMEM164 linked to ferroptosis, autophagy, and Hippo/YAP signalling?
6. Can TMEM164 serve as an early diagnostic biomarker (ROC/LASSO)?

## Repository structure

```
COPD2026/
├── COPD课题/                        # Original CSV data files
├── scripts/
│   ├── 00_setup/
│   │   ├── config.R                 # Global configuration
│   │   └── install_packages.R       # renv-based dependency installation
│   ├── utils/
│   │   ├── theme_nature.R           # Nature publication ggplot2 theme
│   │   ├── geo_helpers.R            # GEO download & parsing
│   │   ├── batch_correction.R       # ComBat-seq, ComBat, Harmony wrappers
│   │   └── marker_lists.R           # Cell-type markers, pathway gene sets
│   ├── 01_data_download/
│   │   ├── download_priority1.R     # 9 core datasets
│   │   ├── download_priority2.R     # 9 validation datasets
│   │   └── download_priority3.R     # 7 extension datasets
│   ├── 02_preprocessing/
│   │   ├── preprocess_bulk_rnaseq.R # DESeq2 + cross-cohort ComBat-seq
│   │   ├── preprocess_microarray.R  # Affy/Agilent/Illumina + ComBat
│   │   ├── preprocess_scrnaseq.R    # Seurat v5 + DoubletFinder + Harmony
│   │   ├── celltype_annotation.R    # AddModuleScore + validation
│   │   └── preprocess_mouse.R       # GSE168299 mouse COPD
│   ├── 03_copd_analysis/
│   │   ├── copd_deg_bulk.R          # Per-cohort DEG + meta-analysis
│   │   ├── copd_deg_sc.R            # Cell-type–resolved DEG
│   │   ├── copd_gold_staging.R      # GOLD 1-4 trajectory + trend test
│   │   ├── copd_marker_interaction.R # Correlation + WGCNA
│   │   └── copd_pathway.R           # GSEA (KEGG/Hallmark) + AUCell
│   ├── 04_ipf_analysis/
│   │   ├── ipf_deg_bulk.R           # IPF DEG meta-analysis
│   │   ├── ipf_deg_sc.R             # KRT8+ transitional + cell-type DEG
│   │   ├── ipf_emt_trajectory.R     # Monocle3/Slingshot AT2→AT1
│   │   ├── ipf_ferroptosis.R        # Ferroptosis/autophagy/Hippo/lipid
│   │   └── ipf_marker_interaction.R # SPP1/MMP7/MUC5B/KRT17/ACTA2
│   ├── 05_integrative/
│   │   ├── cross_disease_comparison.R # COPD vs IPF shared/divergent
│   │   ├── diagnostic_model.R        # ROC + LASSO multi-gene panel
│   │   ├── gene_network.R            # WGCNA + CellChat + TF enrichment
│   │   └── epithelial_fate.R         # 12-axis signalling profiling
│   └── 06_figures/
│       ├── fig1_overview.R           # Pan-cohort + volcano
│       ├── fig2_copd_sc.R            # COPD single-cell UMAP/violin/dot
│       ├── fig3_gold_staging.R       # GOLD staging + ROC
│       ├── fig4_ipf_trajectory.R     # EMT trajectory
│       ├── fig5_pathway.R            # GSEA + AUCell + ferroptosis
│       ├── fig6_cross_disease.R      # Cross-disease integration
│       └── fig_supplementary.R       # QC, batch, mouse, correlation
├── data/                             # (gitignored) raw & processed data
├── results/
│   ├── figures/                      # PDF + PNG outputs
│   ├── tables/                       # CSV result tables
│   └── reports/                      # Session info logs
└── docs/
    └── analysis_report.md            # Full analysis report template
```

## Datasets

### Priority 1 — Core (9 datasets)
| GEO ID | Disease | Type | Platform | Race/Origin |
|--------|---------|------|----------|-------------|
| GSE136831 | COPD+IPF | snRNA-seq | 10x Chromium | Caucasian |
| GSE57148 | COPD | Microarray | Agilent | Korean |
| GSE76925 | COPD | Microarray | Illumina | Caucasian |
| GSE310058 | COPD | Bulk RNA-seq | Illumina HiSeq | Chinese |
| GSE47460 | COPD+IPF | Microarray | Agilent | Caucasian |
| GSE135893 | IPF | scRNA-seq | 10x Chromium | Caucasian |
| GSE134692 | IPF | Bulk RNA-seq | Illumina HiSeq | Japanese |
| GSE124685 | IPF | Bulk RNA-seq | Illumina HiSeq | Australian |
| GSE32537 | IPF | Microarray | Affy HG-U133+ | Caucasian |

### Priority 2 — Validation (9 datasets)
GSE151052, GSE38974, GSE69818, GSE128033, GSE157671, GSE24206, GSE110147, GSE122960, GSE158127

### Priority 3 — Extension (7 datasets)
GSE168299 (mouse), GSE185534, GSE173896, GSE196638, GSE132771, GSE175436, GSE190889

## Quick start

```bash
# 1. Clone
git clone https://github.com/huicmc1223/COPD2026.git
cd COPD2026

# 2. Install dependencies (R ≥ 4.3)
Rscript scripts/00_setup/install_packages.R

# 3. Download data
Rscript scripts/01_data_download/download_priority1.R

# 4. Preprocess
Rscript scripts/02_preprocessing/preprocess_bulk_rnaseq.R
Rscript scripts/02_preprocessing/preprocess_microarray.R
Rscript scripts/02_preprocessing/preprocess_scrnaseq.R
Rscript scripts/02_preprocessing/celltype_annotation.R

# 5. Analysis
Rscript scripts/03_copd_analysis/copd_deg_bulk.R
# ... (see full pipeline below)

# 6. Figures
Rscript scripts/06_figures/fig1_overview.R
```

## Execution order

```
Phase 1: Setup & Download
  install_packages.R → download_priority1.R → download_priority2.R → download_priority3.R

Phase 2: Preprocessing
  preprocess_bulk_rnaseq.R → preprocess_microarray.R → preprocess_scrnaseq.R
  → celltype_annotation.R → preprocess_mouse.R

Phase 3: COPD Analysis
  copd_deg_bulk.R → copd_deg_sc.R → copd_gold_staging.R
  → copd_marker_interaction.R → copd_pathway.R

Phase 4: IPF Analysis
  ipf_deg_bulk.R → ipf_deg_sc.R → ipf_emt_trajectory.R
  → ipf_ferroptosis.R → ipf_marker_interaction.R

Phase 5: Integration
  cross_disease_comparison.R → diagnostic_model.R
  → gene_network.R → epithelial_fate.R

Phase 6: Figures
  fig1_overview.R → fig2_copd_sc.R → fig3_gold_staging.R
  → fig4_ipf_trajectory.R → fig5_pathway.R → fig6_cross_disease.R
  → fig_supplementary.R
```

## Batch correction strategy

| Data type | Method | Package |
|-----------|--------|---------|
| Bulk RNA-seq counts | ComBat-seq | sva |
| Normalized arrays | ComBat | sva |
| scRNA-seq | Harmony | harmony |
| Meta-analysis | REML random effects | metafor |

## Key tools

- **Seurat v5** — scRNA-seq analysis, DoubletFinder, BPCells
- **DESeq2** — Bulk RNA-seq differential expression
- **limma** — Microarray analysis
- **metafor** — Random-effects meta-analysis
- **Monocle3 / Slingshot** — Pseudotime trajectory
- **WGCNA** — Weighted gene co-expression network
- **CellChat** — Cell-cell communication
- **AUCell** — Cell-type–specific pathway activity
- **clusterProfiler** — GSEA (KEGG + MSigDB Hallmark)
- **pROC + glmnet** — Diagnostic biomarker modelling

## Reproducibility

All scripts set `set.seed(2026)` and save `sessionInfo()` to `results/reports/`. Package versions managed via `renv`.

## License

MIT