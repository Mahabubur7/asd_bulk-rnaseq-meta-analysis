# ==============================================================================
# 00_Setup.R — PACKAGE INSTALLATION
# ==============================================================================

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

CRAN_REQUIRED <- c(
  "tidyverse", "rio", "ggplot2", "ggrepel", "ggpubr", "pheatmap", "patchwork",
  "igraph", "ggraph", "scales", "httr", "jsonlite", "Matrix", "metafor",
  "survival", "survminer", "glmnet", "randomForest", "pROC", "caret",
  "VennDiagram", "ggVennDiagram", "msigdbr", "gridExtra",
  "circlize", "future"
)
pak::pkg_install(CRAN_REQUIRED)

# Bioconductor Packages
BIOC_REQUIRED <- c(
  "DESeq2", "sva", "edgeR", "limma", "RUVSeq", "clusterProfiler", "enrichplot",
  "org.Hs.eg.db", "ReactomePA", "GSVA", "recount3", "SummarizedExperiment",
  "MultiAssayExperiment", "curatedTCGAData", "ComplexHeatmap",
  "genekitr", "TCGAbiolinks", "biomaRt", "GEOquery", "ExperimentHub"
)
pak::pkg_install(BIOC_REQUIRED)

# MetaVolcanoR — install from Bioconductor (more reliable than GitHub fork)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("MetaVolcanoR", update = FALSE, ask = FALSE)

# Optional Bioconductor Packages
BIOC_OPTIONAL <- c(
  "decoupleR", "multiMiR", "depmap", "cBioPortalData",
  "SingleR", "celldex", "slingshot", "AUCell"
)
pak::pkg_install(BIOC_OPTIONAL)
