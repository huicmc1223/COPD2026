###############################################################################
# marker_lists.R — Cell-type markers and pathway gene sets
# Derived from:
#   COPD课题/06_Common cell lineages.csv
#   COPD课题/07_Pulmonary Immune cell lineages.csv
#   COPD课题/08_KEGG.csv
# Usage: source(here::here("scripts/utils/marker_lists.R"))
###############################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 1. Common pulmonary cell lineages (from file 06)
# ══════════════════════════════════════════════════════════════════════════════

markers_common <- list(
  AT2 = c("SFTPC", "SFTPB", "NAPSA", "LPCAT1", "SLC34A2", "ABCA3"),
  AT1 = c("AGER", "HOPX", "AQP5", "PDPN", "CAV1", "CLIC5"),
  Basal = c("TP63", "KRT5", "KRT14", "KRT15", "ITGB6"),
  Club = c("SCGB1A1", "SCGB3A2", "SCGB3A1"),
  Goblet = c("MUC5AC", "MUC5B", "SPDEF", "PIGR", "AQP3"),
  Ciliated = c("FOXJ1", "DNAH5", "TPPP3"),
  Fibroblast = c("COL1A1", "COL3A1", "DCN", "PDGFRA", "VIM", "CD34"),
  Myofibroblast = c("ACTA2", "MYLK", "MYH11", "ITGA8", "PDGFRB"),
  AM = c("MARCO", "MRC1", "CD68", "CD163"),
  Endothelial = c("PRX", "FCN3", "CA4", "VWF", "PECAM1"),
  PNEC = c("ASCL1", "CHGA", "CHGB", "CALCA", "GRP", "INSM1", "NEUROD1"),
  Ionocyte = c("FOXI1", "CFTR", "ASCL3", "BSND", "ATP6V1B1", "CLCNKA"),
  Tuft = c("POU2F3", "GNAT3", "TRPM5", "IL25", "DCLK1")
)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Pulmonary immune cell lineages (from file 07)
# ══════════════════════════════════════════════════════════════════════════════

markers_immune <- list(
  AM = c("MARCO", "MRC1", "CD68", "CD163", "FABP4", "PPARG"),
  IM = c("CD163", "MRC1", "FPR3", "LYVE1"),
  Monocyte = c("CD14", "S100A8", "S100A9", "VCAN", "CCR2"),
  DC = c("CLEC9A", "CLEC10A", "CD1C", "CD207"),
  CD4_T = c("CD3D", "CD4", "IL7R", "FOXP3"),
  CD8_T = c("CD3D", "CD8A", "CD8B", "GZMA", "PRF1"),
  B_cell = c("MS4A1", "CD79A"),
  Plasma = c("IGHG4", "IGHG1", "CD79A"),
  NK = c("NCAM1", "KLRD1", "GNLY", "PRF1", "GZMB"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  Neutrophil = c("S100A8", "S100A9", "FCGR3B", "CSF3R")
)

# Combined marker list for annotation
markers_all <- c(markers_common, markers_immune)
# Remove duplicates (AM is in both)
markers_all <- markers_all[!duplicated(names(markers_all))]

# ══════════════════════════════════════════════════════════════════════════════
# 3. Disease-specific markers (from files 06/07)
# ══════════════════════════════════════════════════════════════════════════════

# IPF-specific alterations
markers_ipf_specific <- list(
  IPF_AT2 = c("MUC5B", "SOX2"),
  Atypical_Basal = c("KRT17", "VIM", "CDH2", "FN1", "COL1A1", "TNC",
                      "MMP7", "HMGA2"),
  Profibrotic_AM = c("SPP1", "LPL", "MMP9", "CTSK", "SPARC", "CSF1",
                      "APOE", "TREM2", "CD9", "CD36"),
  KRT8_transitional = c("KRT8", "CLDN4", "LGALS3"),
  EMT_markers = c("VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "TWIST1")
)

# COPD-specific alterations
markers_copd_specific <- list(
  AT2_B_subset = c("HHIP"),
  MT_high_AM = c("MT1G", "MT1X", "MT2A", "HMOX1", "MT1E"),
  Oxidative_stress = c("SOD3", "THBS1"),
  PNEC_COPD = c("OR2W1")
)

# ══════════════════════════════════════════════════════════════════════════════
# 4. Known COPD/IPF biomarkers (from analysis requirements)
# ══════════════════════════════════════════════════════════════════════════════

copd_known_markers <- c("MUC5B", "MMP7", "SPP1", "HHIP", "CC16")
ipf_known_markers  <- c("SPP1", "MMP7", "MUC5B", "KRT17", "ACTA2")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Pathway gene sets (from file 08 — KEGG + epithelial fate)
# ══════════════════════════════════════════════════════════════════════════════

kegg_focus <- list(
  ECM_receptor       = "hsa04512",
  Focal_adhesion     = "hsa04510",
  PI3K_Akt           = "hsa04151",
  TGFb               = "hsa04350",
  Protein_digestion   = "hsa04974",
  MAPK               = "hsa04010",
  Ferroptosis         = "hsa04216",
  Cytokine_receptor   = "hsa04060",
  Chemokine           = "hsa04062",
  NFkB               = "hsa04064",
  TLR                = "hsa04620"
)

# Ferroptosis core genes
ferroptosis_genes <- c("GPX4", "ACSL4", "TFRC", "SLC7A11", "FSP1",
                        "HMOX1", "FTH1", "FTL", "NCOA4", "LPCAT3")

# Autophagy core genes
autophagy_genes <- c("ATG5", "ATG7", "ATG12", "BECN1", "LAMP1", "LAMP2",
                      "MAP1LC3B", "SQSTM1", "ULK1", "MTOR")

# Hippo / YAP pathway
hippo_yap_genes <- c("YAP1", "WWTR1", "TEAD1", "TEAD4", "LATS1", "LATS2",
                      "MST1", "MST2", "MOB1A")

# Epithelial fate reprogramming axes (from file 08 lower table)
epithelial_fate <- list(
  TGFb_axis = c("TGFB1", "SMAD2", "SMAD3", "CTGF", "SERPINE1"),
  Wnt_axis  = c("CTNNB1", "AXIN2", "LGR5"),
  Notch_axis = c("NOTCH1", "JAG1", "HES1"),
  Hippo_axis = c("YAP1", "WWTR1", "TEAD1"),
  p53_axis   = c("TP53", "CDKN1A"),
  ER_stress  = c("ATF4", "DDIT3", "XBP1"),  # DDIT3 = CHOP
  ROS_axis   = c("NFE2L2", "KEAP1"),  # NFE2L2 = NRF2
  NFkB_axis  = c("RELA", "NFKB1"),
  IL6_STAT3  = c("IL6", "STAT3"),
  Glycolysis = c("HK2", "LDHA", "PKM"),
  Lipid_meta = c("SREBF1", "PPARG"),  # TMEM164-related
  Mito_OXPHOS = c("PPARGC1A", "TFAM")
)

# ══════════════════════════════════════════════════════════════════════════════
# 6. Convenience: flatten all markers for quick gene checks
# ══════════════════════════════════════════════════════════════════════════════

all_marker_genes <- unique(unlist(c(
  markers_all, markers_ipf_specific, markers_copd_specific,
  ferroptosis_genes, autophagy_genes, hippo_yap_genes,
  epithelial_fate, "TMEM164"
)))

message("✔ marker_lists.R loaded — ", length(all_marker_genes),
        " unique marker genes registered")
