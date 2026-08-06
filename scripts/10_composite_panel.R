# ==============================================================================
# 10_composite_panel.R — Publication-ready multi-panel figure
# Combines: Heatmap (mean log2FC) + GO-BP GSEA dotplot + PPI network
# ==============================================================================

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(rio)
library(grid)
library(png)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; desc <- dplyr::desc

dir.create("results/figures/composite", showWarnings = FALSE, recursive = TRUE)

# ---- Panel A: Mean log2FC Heatmap ------------------------------------------
meta_degs <- import("results/tables/meta-analysis/filtered_meta_degs.csv")

# Build per-dataset log2FC lookup
deg_files <- list.files("results/tables/DESeq2", pattern = "\\.csv$", full.names = TRUE)
per_dataset <- lapply(deg_files, function(f) {
  df <- import(f)
  df$Gene_ID <- as.character(df$Gene_ID)
  df
})
names(per_dataset) <- tools::file_path_sans_ext(basename(deg_files))

# Get log2FC for each sig gene in each dataset
sig_ids <- as.character(meta_degs$Gene_ID)
lfc_matrix <- sapply(per_dataset, function(df) {
  matched <- df$log2FoldChange[match(sig_ids, df$Gene_ID)]
  matched
})
rownames(lfc_matrix) <- sig_ids

# Add gene symbols as row labels
gene_labels <- sig_ids
anno_files <- list.files("results/tables/annotated", pattern = "\\.csv$", full.names = TRUE)
if (length(anno_files) > 0) {
  gene_anno <- lapply(anno_files, function(f) {
    import(f) |> select(any_of(c("Gene_ID", "Gene_Symbol")))
  }) |> bind_rows() |> distinct(Gene_ID, .keep_all = TRUE)
  matched_sym <- gene_anno$Gene_Symbol[match(sig_ids, gene_anno$Gene_ID)]
  gene_labels <- ifelse(!is.na(matched_sym) & matched_sym != "", matched_sym, sig_ids)
}
rownames(lfc_matrix) <- gene_labels

# Mean log2FC and direction consistency
mean_lfc <- rowMeans(lfc_matrix, na.rm = TRUE)
direction_consistent <- apply(lfc_matrix, 1, function(row) {
  row <- row[!is.na(row)]
  if (length(row) < 2) return("Consistent")
  if (all(row > 0) || all(row < 0)) "Consistent" else "Inconsistent"
})
meta_status <- ifelse(mean_lfc > 0, "Up-regulated", "Down-regulated")

# Row annotations
row_anno <- rowAnnotation(
  Direction_Consistent = direction_consistent,
  Meta_Status = meta_status,
  col = list(
    Direction_Consistent = c("Consistent" = "#1b9e77", "Inconsistent" = "#d95f02"),
    Meta_Status = c("Up-regulated" = "#e34a33", "Down-regulated" = "#2c7fb8")
  ),
  annotation_name_gp = gpar(fontsize = 8)
)

col_fun_lfc <- colorRamp2(c(-1.5, 0, 1.5), c("#2166ac", "white", "#b2182b"))

ht <- Heatmap(
  as.matrix(mean_lfc),
  name = "Mean\nlog2FC",
  col = col_fun_lfc,
  right_annotation = row_anno,
  row_labels = gene_labels,
  row_names_gp = gpar(fontsize = 7),
  show_column_names = TRUE,
  column_labels = "Meta_Mean_log2FC",
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  width = unit(1.5, "cm"),
  column_title = "Meta-Analysis DEGs (Mean log2FC)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(title = "Mean\nlog2FC", direction = "vertical",
                               title_gp = gpar(fontsize = 8),
                               labels_gp = gpar(fontsize = 7))
)

# Save heatmap as standalone PNG for compositing
png("results/figures/composite/_panel_A_heatmap.png",
    width = 5, height = 10, units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# Also save mean log2FC heatmap to heatmap folder
png("results/figures/heatmap/meta_degs_mean_log2FC_heatmap.png",
    width = 5, height = 10, units = "in", res = 300)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf("results/figures/heatmap/meta_degs_mean_log2FC_heatmap.pdf",
    width = 5, height = 10)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Panel A (heatmap) saved")

# ---- Panel B: GO-BP GSEA dotplot -------------------------------------------
gobp_file <- "results/tables/enrichment/asd_ipsc_gsea_GSEA_GO_BP.csv"
panel_B <- NULL
if (file.exists(gobp_file)) {
  gobp <- import(gobp_file)
  gobp <- gobp[order(gobp$p.adjust), ]
  gobp <- head(gobp, 15)
  gobp$Label <- stringr::str_wrap(gsub("_", " ", tools::toTitleCase(tolower(gobp$Description))), 30)

  panel_B <- ggplot(gobp, aes(x = NES, y = reorder(Label, NES))) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES, yend = reorder(Label, NES)),
                 color = "grey80", linewidth = 0.5) +
    geom_point(aes(size = setSize, fill = p.adjust), shape = 21, color = "grey30") +
    scale_fill_gradient(low = "#b2182b", high = "#2166ac", trans = "log10", name = "p.adjust") +
    scale_size_continuous(range = c(3, 8), name = "Set size") +
    labs(title = "GSEA: GO Biological Process", x = "NES", y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5))

  ggsave("results/figures/composite/_panel_B_gobp.png",
         plot = panel_B, width = 7, height = 6, dpi = 300, bg = "white")
  message("Panel B (GO-BP) saved")
} else {
  message("GO-BP GSEA file not found, skipping panel B")
}

# ---- Panel C: PPI network (use existing ggraph) ---------------------------
ppi_file <- "results/figures/network/asd_ipsc_ppi_ggraph.png"
panel_C_exists <- file.exists(ppi_file)
if (panel_C_exists) {
  message("Panel C (PPI network) exists at ", ppi_file)
} else {
  message("PPI ggraph not found, skipping panel C")
}

# ---- Composite: combine all three panels -----------------------------------
# Read panels as raster images and arrange
library(gridExtra)

img_A <- readPNG("results/figures/composite/_panel_A_heatmap.png")
grob_A <- rasterGrob(img_A, interpolate = TRUE)

panels <- list(grob_A)
widths <- c(2)
labels <- c("A")

if (!is.null(panel_B)) {
  ggsave("results/figures/composite/_panel_B_gobp.png",
         plot = panel_B, width = 7, height = 10, dpi = 300, bg = "white")
  img_B <- readPNG("results/figures/composite/_panel_B_gobp.png")
  grob_B <- rasterGrob(img_B, interpolate = TRUE)
  panels <- c(panels, list(grob_B))
  widths <- c(widths, 3)
  labels <- c(labels, "B")
}

if (panel_C_exists) {
  img_C <- readPNG(ppi_file)
  grob_C <- rasterGrob(img_C, interpolate = TRUE)
  panels <- c(panels, list(grob_C))
  widths <- c(widths, 3)
  labels <- c(labels, "C")
}

# Add panel labels
labeled_panels <- lapply(seq_along(panels), function(i) {
  arrangeGrob(
    textGrob(labels[i], x = 0.05, y = 0.95, just = c("left", "top"),
             gp = gpar(fontsize = 18, fontface = "bold")),
    panels[[i]],
    ncol = 1, heights = c(0.05, 0.95)
  )
})

composite <- arrangeGrob(grobs = labeled_panels, ncol = length(panels),
                          widths = widths)

png("results/figures/composite/Figure1_Composite_Panel.png",
    width = 20, height = 12, units = "in", res = 300)
grid.draw(composite)
dev.off()

pdf("results/figures/composite/Figure1_Composite_Panel.pdf",
    width = 20, height = 12)
grid.draw(composite)
dev.off()

message("\nComposite figure saved to results/figures/composite/ (PNG + PDF)")
message("  Panels: ", paste(labels, collapse = " + "))
