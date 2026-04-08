# Analysis Report: TMEM164 in COPD & IPF

> **Multi-cohort transcriptomic profiling of TMEM164 as a COPD early diagnostic biomarker and IPF fibrosis mediator**

---

## 1. Study design

### 1.1 Rationale
TMEM164 (Transmembrane protein 164) has emerging evidence linking it to lipid metabolism, ferroptosis, and epithelial cell fate — processes central to both COPD pathogenesis and IPF fibrosis. This study leverages 25 public GEO datasets spanning diverse platforms and ethnicities to establish TMEM164 as a robust biomarker candidate.

### 1.2 Datasets
- **18 COPD datasets**: bulk RNA-seq, microarray (Affy/Agilent/Illumina), scRNA-seq, snRNA-seq
- **16 IPF datasets**: overlapping platforms with COPD, plus spatial transcriptomics
- **1 mouse model**: GSE168299 (Tmem164 ortholog)
- **Ethnicities**: Korean, Caucasian, Japanese, Chinese, Australian

### 1.3 Analysis pipeline
Executed in 6 phases: Setup → Download → Preprocessing → Disease-specific → Integration → Figures.

---

## 2. Preprocessing

### 2.1 Bulk RNA-seq
- DESeq2 pipeline: size-factor normalization → variance-stabilising transform (VST)
- QC: PCA + sample correlation heatmap
- Cross-cohort batch correction: ComBat-seq (within-disease)

### 2.2 Microarray
- Platform-specific: Affy (RMA via oligo), Agilent (quantile via limma), Illumina (neqc via limma)
- Probe-to-gene: max-mean collapse
- Cross-platform: ComBat on quantile-normalized expression

### 2.3 scRNA-seq / snRNA-seq
- Per-dataset Seurat v5: SCTransform v2, DoubletFinder (expected 5%)
- Multi-dataset integration: Harmony (θ = 2, max_iter = 30)
- BPCells/HDF5 for memory management of large datasets
- QC thresholds: 200 ≤ nFeature ≤ 8000, percent.mt ≤ 20%

### 2.4 Cell-type annotation
- Method: AddModuleScore (13 common + 11 immune + disease-specific lineages)
- Validation: DotPlot of canonical markers, UMAP overlay, per-cluster composition

---

## 3. COPD: TMEM164 differential expression

### 3.1 Bulk / microarray meta-analysis
- Per-cohort DEG: DESeq2 (RNA-seq) or limma (array)
- Meta-analysis: metafor::rma() with REML estimation
- **Result**: Pooled log2FC = _[to be filled]_, 95% CI = _[to be filled]_, I² = _[to be filled]_
- Forest plot: `results/figures/copd_tmem164_forest.pdf`

### 3.2 Single-cell resolution
- Per-cell-type FindMarkers (Wilcoxon, COPD vs Control)
- Cell types showing significant TMEM164 change: _[to be filled]_
- Visualisation: violin plot, DotPlot, FeaturePlot

### 3.3 GOLD staging trajectory
- Datasets: GSE310058, GSE38974, GSE76925, GSE57148
- Jonckheere-Terpstra trend test: p = _[to be filled]_
- Direction: _[monotonic increase / decrease / non-linear]_

---

## 4. IPF: TMEM164 and fibrosis

### 4.1 Bulk DEG meta-analysis
- Pooled log2FC = _[to be filled]_
- Forest plot: `results/figures/ipf_tmem164_forest.pdf`

### 4.2 KRT8+ transitional state
- Identification: top 10% KRT8 transitional module score
- TMEM164 enrichment in KRT8+ vs other epithelial: log2FC = _[to be filled]_

### 4.3 EMT trajectory
- Method: Monocle3 (primary) / Slingshot (fallback)
- Root: AT2 cells, trajectory: AT2 → KRT8+ → AT1
- TMEM164 dynamics: _[transient peak / progressive increase / etc.]_
- IPF vs Control divergence: `results/figures/ipf_emt_trajectory_tmem164.pdf`

### 4.4 Ferroptosis / autophagy / Hippo
- Key correlations with TMEM164:
  - GPX4 (ρ = _[to be filled]_)
  - ACSL4 (ρ = _[to be filled]_)
  - BECN1 (ρ = _[to be filled]_)
  - YAP1 (ρ = _[to be filled]_)

---

## 5. Integrative analysis

### 5.1 Cross-disease comparison
- Shared regulation: cell types where TMEM164 is dysregulated in both COPD and IPF
- Divergent patterns: cell types with opposing direction
- Figure: `results/figures/cross_disease_deg_comparison.pdf`

### 5.2 Diagnostic model
- Single-gene ROC (COPD vs Control): AUC = _[to be filled]_ per cohort
- LASSO multi-gene panel: AUC = _[to be filled]_
- Selected features: _[to be filled]_

### 5.3 Gene regulatory network
- WGCNA module containing TMEM164: module _[to be filled]_ (_[N]_ genes)
- Top hub genes: _[to be filled]_
- GO enrichment of hub genes: top terms = _[to be filled]_

### 5.4 Epithelial fate axes
- Significant pathways (TMEM164-high vs -low):
  - COPD: _[to be filled]_
  - IPF: _[to be filled]_

---

## 6. Mouse validation

- GSE168299: Tmem164 in COPD mouse model
- Direction consistent with human: _[yes/no]_
- Figure: `results/figures/supp_S6_mouse_tmem164.pdf`

---

## 7. Summary of key findings

| Finding | Evidence | Strength |
|---------|----------|----------|
| TMEM164 is differentially expressed in COPD | Meta-analysis across N cohorts | _[to be filled]_ |
| Expression correlates with GOLD severity | Trend test | _[to be filled]_ |
| AT2 cells show strongest change | scRNA-seq cell-type DEG | _[to be filled]_ |
| TMEM164 marks EMT transition in IPF | Pseudotime analysis | _[to be filled]_ |
| Association with ferroptosis | Spearman correlation | _[to be filled]_ |
| Diagnostic potential (AUC > 0.7) | ROC analysis | _[to be filled]_ |

---

## 8. Reproducibility

- All random seeds: `2026`
- R version: _[to be filled]_
- Session info: `results/reports/`
- Package management: `renv`

---

*Report generated by the COPD2026 analysis pipeline. Fill placeholders after running the complete analysis.*
