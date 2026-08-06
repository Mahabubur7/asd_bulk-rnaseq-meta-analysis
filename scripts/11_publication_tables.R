# ==============================================================================
# 11_publication_tables.R — Publication-ready tables for manuscript
# Generates formatted tables using gt (display) and flextable (Word export)
# ==============================================================================

# Install if needed
if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")

library(tidyverse)
library(rio)
library(gt)
library(flextable)
library(officer)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; arrange <- dplyr::arrange; desc <- dplyr::desc

dir.create("results/tables/publication", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/tables", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# TABLE 1: Dataset Characteristics
# ==============================================================================

table1_df <- tibble(
  Attribute = c("GEO Accession", "Source Study", "Model System", "Comparison",
                "Total Samples", "ASD Samples", "Control Samples", "Platform",
                "SVA Surrogate Variables", "Protein-Coding Genes Retained"),
  GSE67528 = c("GSE67528", "Marchetto et al., 2017", 
               "iPSC + iPSC-derived NPCs + neurons", "ASD vs. Control",
               "83", "55", "28", "RNA-Seq (Illumina)",
               "8", "16,852"),
  GSE124308 = c("GSE124308", "DeRosa et al., 2018",
                "iPSC-derived cortical neurons", "ASD vs. Control",
                "22", "12", "10", "RNA-Seq (Illumina)",
                "4", "16,365")
)

# gt version (for display / HTML / PDF)
gt_table1 <- table1_df |>
  gt() |>
  tab_header(
    title = md("**Table 1.** Dataset characteristics"),
    subtitle = "Two independent iPSC-based ASD RNA-Seq datasets included in the meta-analysis"
  ) |>
  cols_label(
    Attribute = "Attribute",
    GSE67528 = md("**GSE67528**"),
    GSE124308 = md("**GSE124308**")
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Attribute)
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#1F3864"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_labels()
  ) |>
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(11),
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    column_labels.border.bottom.style = "solid"
  )

gtsave(gt_table1, "results/figures/tables/Table1_dataset_characteristics.png", vwidth = 900)
gtsave(gt_table1, "results/figures/tables/Table1_dataset_characteristics.html")

# flextable version (for Word export)
ft_table1 <- flextable(table1_df) |>
  set_header_labels(Attribute = "Attribute", GSE67528 = "GSE67528", GSE124308 = "GSE124308") |>
  bold(j = 1) |>
  bg(part = "header", bg = "#1F3864") |>
  color(part = "header", color = "white") |>
  bold(part = "header") |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Calibri", part = "all") |>
  autofit() |>
  border_outer(border = fp_border(color = "black", width = 1.5)) |>
  border_inner_h(border = fp_border(color = "grey80", width = 0.5)) |>
  set_caption("Table 1. Characteristics of the two datasets included in the meta-analysis.")

save_as_docx(ft_table1, path = "results/tables/publication/Table1_dataset_characteristics.docx")
message("Table 1 saved")


# ==============================================================================
# TABLE 2: Meta-Analysis DEGs (full list with regulation)
# ==============================================================================

meta_degs <- import("results/tables/meta-analysis/filtered_meta_degs.csv")

table2_df <- meta_degs |>
  filter(!is.na(Gene_Symbol), Gene_Symbol != "") |>
  mutate(
    Regulation = ifelse(randomSummary > 0, "Up", "Down"),
    `log2FC` = round(randomSummary, 2),
    `P-value` = formatC(randomP, format = "e", digits = 2),
    `CI Lower` = round(randomCi.lb, 2),
    `CI Upper` = round(randomCi.ub, 2)
  ) |>
  select(Gene_Symbol, Gene_Description, `log2FC`, `P-value`, `CI Lower`, `CI Upper`, Regulation) |>
  arrange(desc(abs(`log2FC`)))

gt_table2 <- table2_df |>
  gt() |>
  tab_header(
    title = md("**Table 2.** Significant meta-analysis DEGs"),
    subtitle = "Random-effects model: randomP < 0.05, |log2FC| ≥ 1"
  ) |>
  cols_label(
    Gene_Symbol = "Gene",
    Gene_Description = "Description",
    `log2FC` = "log2FC",
    `P-value` = "P-value",
    `CI Lower` = "95% CI Lower",
    `CI Upper` = "95% CI Upper",
    Regulation = "Direction"
  ) |>
  tab_style(
    style = cell_text(color = "#de2d26", weight = "bold"),
    locations = cells_body(columns = Regulation, rows = Regulation == "Up")
  ) |>
  tab_style(
    style = cell_text(color = "#2171b5", weight = "bold"),
    locations = cells_body(columns = Regulation, rows = Regulation == "Down")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#1F3864"), cell_text(color = "white", weight = "bold")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = Gene_Symbol)
  ) |>
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(14),
    data_row.padding = px(4)
  )

gtsave(gt_table2, "results/figures/tables/Table2_meta_degs.png", vwidth = 1100)
gtsave(gt_table2, "results/figures/tables/Table2_meta_degs.html")

ft_table2 <- flextable(table2_df) |>
  set_header_labels(
    Gene_Symbol = "Gene", Gene_Description = "Description",
    `log2FC` = "log2FC", `P-value` = "P-value",
    `CI Lower` = "95% CI Lower", `CI Upper` = "95% CI Upper",
    Regulation = "Direction"
  ) |>
  italic(j = 1) |>
  bg(part = "header", bg = "#1F3864") |>
  color(part = "header", color = "white") |>
  bold(part = "header") |>
  color(i = ~ Regulation == "Up", j = "Regulation", color = "#de2d26") |>
  bold(i = ~ Regulation == "Up", j = "Regulation") |>
  color(i = ~ Regulation == "Down", j = "Regulation", color = "#2171b5") |>
  bold(i = ~ Regulation == "Down", j = "Regulation") |>
  fontsize(size = 9, part = "all") |>
  font(fontname = "Calibri", part = "all") |>
  autofit() |>
  border_outer(border = fp_border(color = "black", width = 1.5)) |>
  border_inner_h(border = fp_border(color = "grey80", width = 0.5)) |>
  set_caption("Table 2. Significant differentially expressed genes identified by random-effects meta-analysis (randomP < 0.05, |log2FC| ≥ 1). Genes are listed by decreasing absolute fold-change.")

save_as_docx(ft_table2, path = "results/tables/publication/Table2_meta_degs.docx")
message("Table 2 saved")


# ==============================================================================
# TABLE 3: DEG Summary Statistics
# ==============================================================================

table3_df <- tibble(
  Direction = c("Upregulated in ASD", "Downregulated in ASD", "Total"),
  Count = c(18, 21, 39),
  `% of Total Genes` = c("0.089%", "0.103%", "0.192%"),
  `Top Genes` = c(
    "F8A3, DDR1, TUBB, RING1, RXRB, HLA-A",
    "GABBR1, DDAH2, COL11A2, SP140L, CSNK2B, NNAT",
    ""
  )
)

gt_table3 <- table3_df |>
  gt() |>
  tab_header(
    title = md("**Table 3.** Summary of meta-analysis DEGs by regulation direction")
  ) |>
  tab_style(
    style = list(cell_fill(color = "#1F3864"), cell_text(color = "white", weight = "bold")),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Direction == "Total")
  ) |>
  tab_options(table.font.size = px(12), heading.title.font.size = px(14))

gtsave(gt_table3, "results/figures/tables/Table3_deg_summary.png", vwidth = 900)

ft_table3 <- flextable(table3_df) |>
  bg(part = "header", bg = "#1F3864") |>
  color(part = "header", color = "white") |>
  bold(part = "header") |>
  bold(i = 3) |>
  fontsize(size = 10, part = "all") |>
  font(fontname = "Calibri", part = "all") |>
  autofit() |>
  border_outer(border = fp_border(color = "black", width = 1.5)) |>
  border_inner_h(border = fp_border(color = "grey80", width = 0.5)) |>
  set_caption("Table 3. Summary of meta-analysis DEGs by regulation direction.")

save_as_docx(ft_table3, path = "results/tables/publication/Table3_deg_summary.docx")
message("Table 3 saved")


# ==============================================================================
# TABLE 4: Hub Genes from PPI Network
# ==============================================================================

hub_file <- "results/tables/network/asd_ipsc_hub_genes.csv"
if (file.exists(hub_file)) {
  hub_raw <- import(hub_file)
  
  table4_df <- hub_raw |>
    mutate(
      Degree = as.integer(Degree),
      Betweenness = round(Betweenness, 3),
      Closeness = round(Closeness, 3),
      Eigenvector = round(Eigenvector, 3),
      `log2FC` = round(Log2FC, 2),
      Direction = ifelse(Log2FC > 0, "Up", "Down"),
      Function = case_when(
        Symbol == "HLA-A" ~ "MHC class I; synaptic pruning, neuroinflammation",
        Symbol == "RING1" ~ "Polycomb PRC1; chromatin remodeling",
        Symbol == "PRRC2A" ~ "RNA-binding; m6A reader, brain-enriched",
        Symbol == "TRIM26" ~ "E3 ubiquitin ligase; MHC locus, immune signaling",
        TRUE ~ ""
      )
    ) |>
    select(Symbol, Degree, Betweenness, Closeness, Eigenvector, `log2FC`, Direction, Function)
  
  gt_table4 <- table4_df |>
    gt() |>
    tab_header(
      title = md("**Table 4.** Hub genes identified by PPI network centrality analysis"),
      subtitle = "Nodes with degree centrality in the top 10th percentile"
    ) |>
    cols_label(Symbol = "Gene") |>
    tab_style(
      style = list(cell_fill(color = "#1F3864"), cell_text(color = "white", weight = "bold")),
      locations = cells_column_labels()
    ) |>
    tab_style(
      style = cell_text(style = "italic", weight = "bold"),
      locations = cells_body(columns = Symbol)
    ) |>
    tab_style(
      style = cell_text(color = "#de2d26", weight = "bold"),
      locations = cells_body(columns = Direction)
    ) |>
    tab_options(table.font.size = px(11), heading.title.font.size = px(14))
  
  gtsave(gt_table4, "results/figures/tables/Table4_hub_genes.png", vwidth = 1100)
  
  ft_table4 <- flextable(table4_df) |>
    set_header_labels(Symbol = "Gene") |>
    italic(j = 1) |> bold(j = 1) |>
    bg(part = "header", bg = "#1F3864") |>
    color(part = "header", color = "white") |>
    bold(part = "header") |>
    color(j = "Direction", color = "#de2d26") |>
    bold(j = "Direction") |>
    fontsize(size = 10, part = "all") |>
    font(fontname = "Calibri", part = "all") |>
    autofit() |>
    border_outer(border = fp_border(color = "black", width = 1.5)) |>
    border_inner_h(border = fp_border(color = "grey80", width = 0.5)) |>
    set_caption("Table 4. Hub genes identified by PPI network centrality analysis (top 10% degree).")
  
  save_as_docx(ft_table4, path = "results/tables/publication/Table4_hub_genes.docx")
  message("Table 4 saved")
}


# ==============================================================================
# TABLE 5: Network Modules
# ==============================================================================

modules_file <- "results/tables/network/asd_ipsc_network_modules.csv"
if (file.exists(modules_file)) {
  mod_raw <- import(modules_file)
  
  table5_df <- mod_raw |>
    group_by(Module) |>
    summarise(
      Genes = paste(Symbol, collapse = ", "),
      Size = n(),
      `Mean log2FC` = round(mean(Log2FC), 2),
      .groups = "drop"
    ) |>
    mutate(
      `Functional Theme` = case_when(
        Module == 1 ~ "MHC / immune – cytoskeletal",
        Module == 2 ~ "Chromatin – nuclear receptor signaling",
        Module == 3 ~ "RNA processing – immune",
        TRUE ~ ""
      )
    ) |>
    select(Module, Genes, Size, `Mean log2FC`, `Functional Theme`)
  
  gt_table5 <- table5_df |>
    gt() |>
    tab_header(
      title = md("**Table 5.** PPI network community modules (Louvain algorithm)")
    ) |>
    tab_style(
      style = list(cell_fill(color = "#1F3864"), cell_text(color = "white", weight = "bold")),
      locations = cells_column_labels()
    ) |>
    tab_style(
      style = cell_text(style = "italic"),
      locations = cells_body(columns = Genes)
    ) |>
    tab_options(table.font.size = px(11), heading.title.font.size = px(14))
  
  gtsave(gt_table5, "results/figures/tables/Table5_network_modules.png", vwidth = 1000)
  
  ft_table5 <- flextable(table5_df) |>
    italic(j = "Genes") |>
    bg(part = "header", bg = "#1F3864") |>
    color(part = "header", color = "white") |>
    bold(part = "header") |>
    fontsize(size = 10, part = "all") |>
    font(fontname = "Calibri", part = "all") |>
    autofit() |>
    border_outer(border = fp_border(color = "black", width = 1.5)) |>
    border_inner_h(border = fp_border(color = "grey80", width = 0.5)) |>
    set_caption("Table 5. PPI network community modules detected by the Louvain algorithm.")
  
  save_as_docx(ft_table5, path = "results/tables/publication/Table5_network_modules.docx")
  message("Table 5 saved")
}


# ==============================================================================
# COMBINED: All tables in one Word document
# ==============================================================================

all_tables_doc <- read_docx()

all_tables_doc <- all_tables_doc |>
  body_add_par("Publication Tables", style = "heading 1") |>
  body_add_par("") |>
  body_add_par("Table 1. Characteristics of the two datasets included in the meta-analysis.", style = "heading 2") |>
  body_add_flextable(ft_table1) |>
  body_add_break() |>
  body_add_par("Table 2. Significant differentially expressed genes identified by random-effects meta-analysis.", style = "heading 2") |>
  body_add_flextable(ft_table2) |>
  body_add_break() |>
  body_add_par("Table 3. Summary of meta-analysis DEGs by regulation direction.", style = "heading 2") |>
  body_add_flextable(ft_table3) |>
  body_add_break() |>
  body_add_par("Table 4. Hub genes identified by PPI network centrality analysis.", style = "heading 2")

if (exists("ft_table4")) {
  all_tables_doc <- all_tables_doc |> body_add_flextable(ft_table4) |> body_add_break()
}

all_tables_doc <- all_tables_doc |>
  body_add_par("Table 5. PPI network community modules.", style = "heading 2")

if (exists("ft_table5")) {
  all_tables_doc <- all_tables_doc |> body_add_flextable(ft_table5)
}

print(all_tables_doc, target = "results/tables/publication/All_Publication_Tables.docx")

message("\n=== All publication tables generated ===")
message("Individual tables: results/tables/publication/Table1-5_*.docx")
message("Combined document: results/tables/publication/All_Publication_Tables.docx")
message("PNG versions: results/figures/tables/Table1-5_*.png")
