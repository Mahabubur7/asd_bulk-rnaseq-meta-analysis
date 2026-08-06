# ==============================================================================
# 09_heatmap.R — Heatmap of significant meta-analysis DEGs
# Autism Spectrum Disorder iPSC models (ASD vs Healthy Control)
# ==============================================================================

library(ComplexHeatmap)
library(circlize)
library(DESeq2)
library(rio)
library(tidyverse)

# Pin dplyr verbs
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate

# Output directories
dir.create("results/figures/heatmap", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load meta-analysis significant DEGs --------------------------------
meta_degs <- import("results/tables/meta-analysis/filtered_meta_degs.csv")
message("Significant meta-analysis DEGs to plot: ", nrow(meta_degs))

if (nrow(meta_degs) == 0) {
  stop("No significant DEGs found in filtered_meta_degs.csv — nothing to plot.")
}

# ---- 2. Load raw counts + metadata for both datasets -----------------------
geo_ids <- c("GSE67528", "GSE124308")

load_and_normalize <- function(geo_id) {
  count_file <- paste0("data/raw_counts/", geo_id, "_raw_counts.tsv")
  if (!file.exists(count_file)) {
    count_file <- paste0("data/raw_counts/", geo_id, "_raw_counts.csv")
  }
  count_data_raw <- import(count_file)
  colnames(count_data_raw)[1] <- "GeneID"

  metadata <- import(paste0("data/metadata/", geo_id, "_metadata.csv"))

  gene_ids <- count_data_raw$GeneID
  count_data <- count_data_raw |>
    column_to_rownames("GeneID") |>
    select(any_of(metadata$sample))
  # Round first (GSE67528 has RSEM-style non-integer counts), then coerce
  count_data <- as.data.frame(lapply(count_data, function(x) as.integer(round(as.numeric(x)))))
  rownames(count_data) <- gene_ids[seq_len(nrow(count_data))]

  metadata <- metadata |>
    filter(.data$sample %in% colnames(count_data)) |>
    arrange(match(.data$sample, colnames(count_data)))

  # Standardize condition labels
  metadata$condition <- ifelse(
    grepl("control|healthy|normal", tolower(metadata$condition)),
    "Control", "ASD"
  )

  colData <- data.frame(
    condition = factor(metadata$condition, levels = c("Control", "ASD")),
    row.names = colnames(count_data)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData   = colData,
    design    = ~ condition
  )

  # Variance-stabilizing transformation for visualization
  vsd <- vst(dds, blind = TRUE)

  list(
    vsd       = assay(vsd),
    condition = colData$condition,
    geo_id    = geo_id
  )
}

datasets <- lapply(geo_ids, load_and_normalize)
names(datasets) <- geo_ids

# ---- 3. Merge normalized expression across datasets ------------------------
common_genes <- Reduce(base::intersect, lapply(datasets, function(d) rownames(d$vsd)))
sig_genes    <- meta_degs$Gene_ID[meta_degs$Gene_ID %in% common_genes]

if (length(sig_genes) == 0) {
  stop("No significant DEGs found in common across both datasets.")
}

message("Genes in heatmap: ", length(sig_genes))

# Combine expression matrices
combined_expr <- do.call(cbind, lapply(datasets, function(d) d$vsd[sig_genes, ]))

# Z-score scale per gene (row) for visualization
mat <- t(scale(t(combined_expr)))

# ---- 4. Build annotation bars ----------------------------------------------
condition_vec <- unlist(lapply(datasets, function(d) as.character(d$condition)))
dataset_vec <- unlist(lapply(datasets, function(d) rep(d$geo_id, ncol(d$vsd))))

col_anno <- HeatmapAnnotation(
  Condition = condition_vec,
  Dataset   = dataset_vec,
  col = list(
    Condition = c("ASD" = "#de2d26", "Control" = "#2171b5"),
    Dataset   = c("GSE67528" = "#66c2a5", "GSE124308" = "#fc8d62")
  ),
  annotation_name_side = "left"
)

# ---- 5. Add gene symbol labels if available --------------------------------
gene_labels <- sig_genes
anno_file <- list.files("results/tables/annotated", pattern = "\\.csv$", full.names = TRUE)
if (length(anno_file) > 0) {
  gene_anno <- lapply(anno_file, function(f) {
    import(f) |> select(any_of(c("Gene_ID", "Gene_Symbol")))
  }) |>
    bind_rows() |>
    distinct(Gene_ID, .keep_all = TRUE)

  matched <- gene_anno$Gene_Symbol[match(sig_genes, gene_anno$Gene_ID)]
  gene_labels <- ifelse(!is.na(matched) & matched != "", matched, sig_genes)
}

# ---- 6. Draw heatmap -------------------------------------------------------
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))

ht <- Heatmap(
  mat,
  name = "Z-score",
  col  = col_fun,
  top_annotation = col_anno,
  row_labels     = gene_labels,
  row_names_gp   = gpar(fontsize = if (length(sig_genes) > 50) 5 else 8),
  show_column_names = FALSE,
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  show_row_dend   = TRUE,
  show_column_dend = TRUE,
  column_title    = "Meta-Analysis Significant DEGs (ASD vs Control)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title = "Z-score", direction = "vertical")
)

# Save
png("results/figures/heatmap/meta_degs_heatmap.png",
    width = 10, height = max(6, length(sig_genes) * 0.15 + 3),
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Also save as PDF for publication
pdf("results/figures/heatmap/meta_degs_heatmap.pdf",
    width = 10, height = max(6, length(sig_genes) * 0.15 + 3))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Heatmap saved to results/figures/heatmap/ (PNG + PDF)")
message("  Rows: ", length(sig_genes), " significant DEGs")
message("  Columns: ", ncol(mat), " samples across ", length(geo_ids), " datasets")
