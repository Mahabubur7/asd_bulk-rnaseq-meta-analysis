# Functional enrichment (ORA + GSEA) of ASD iPSC meta-analysis DEGs (ASD vs Healthy Control)

library(tidyverse)
library(rio)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)

options(timeout = 300)  # KEGG REST API downloads can exceed the 60s default

# Pin dplyr verbs
select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
rename <- dplyr::rename; arrange <- dplyr::arrange; desc <- dplyr::desc

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

META_FILE <- "results/tables/meta-analysis/random_effect_model.csv"
OUT_CSV   <- "results/tables/enrichment"
OUT_FIG   <- "results/figures/enrichment"
PREFIX    <- "asd_ipsc"
dir.create(OUT_CSV, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

# --- ID conversion helper ---------------------------------------------------
# Gene_IDs in the meta-analysis output are Ensembl (ENSG...).
# enrichGO/gseGO accept keyType = "ENSEMBL" directly.
# enrichKEGG, enrichPathway, gseKEGG, and Hallmark GSEA require Entrez IDs.
ensembl_to_entrez <- function(ensembl_ids) {
  mapped <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  mapped
}

build_ranked_list <- function(meta_df, id_col = "Gene_ID") {
  df <- meta_df |>
    filter(!is.na(randomP), !is.na(randomSummary), randomP > 0, !is.na(.data[[id_col]])) |>
    mutate(rank_score = sign(randomSummary) * -log10(randomP),
           !!id_col := as.character(.data[[id_col]])) |>
    filter(!duplicated(.data[[id_col]])) |>
    arrange(desc(rank_score))
  setNames(df$rank_score, df[[id_col]])
}

run_ora <- function(ensembl_ids, universe_ensembl) {
  if (length(ensembl_ids) < 5) {
    message("    Too few genes for ORA (n=", length(ensembl_ids), "), skipping.")
    return(NULL)
  }
  results <- list()
  
  # GO: use Ensembl IDs directly with keyType = "ENSEMBL"
  for (ont in c("BP", "MF", "CC")) {
    results[[paste0("GO_", ont)]] <- tryCatch(
      enrichGO(gene = ensembl_ids, universe = universe_ensembl, OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL", ont = ont, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) { message("    GO_", ont, " error: ", e$message); NULL }
    )
  }
  
  # KEGG & Reactome: need Entrez IDs
  entrez_map <- tryCatch(ensembl_to_entrez(ensembl_ids), error = function(e) NULL)
  entrez_universe <- tryCatch(ensembl_to_entrez(universe_ensembl), error = function(e) NULL)
  
  if (!is.null(entrez_map) && nrow(entrez_map) >= 5) {
    eid <- unique(entrez_map$ENTREZID)
    uid <- if (!is.null(entrez_universe)) unique(entrez_universe$ENTREZID) else NULL
    
    results[["KEGG"]] <- tryCatch(
      enrichKEGG(gene = eid, universe = uid, organism = "hsa",
                 pAdjustMethod = "BH", pvalueCutoff = 0.05),
      error = function(e) { message("    KEGG error: ", e$message); NULL }
    )
    
    results[["Reactome"]] <- tryCatch(
      enrichPathway(gene = eid, universe = uid, organism = "human",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
      error = function(e) { message("    Reactome error: ", e$message); NULL }
    )
  } else {
    message("    Entrez conversion yielded too few genes; skipping KEGG/Reactome ORA.")
  }
  
  results
}

run_gsea <- function(ranked_ensembl) {
  if (length(ranked_ensembl) < 100) {
    message("    Too few ranked genes (n=", length(ranked_ensembl), "), skipping GSEA.")
    return(NULL)
  }
  results <- list()
  
  # GO-BP GSEA: Ensembl directly
  results[["GSEA_GO_BP"]] <- tryCatch(
    gseGO(geneList = ranked_ensembl, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
          ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
    error = function(e) { message("    GSEA GO-BP error: ", e$message); NULL }
  )
  
  # Convert ranked list to Entrez for KEGG and Hallmark
  entrez_map <- tryCatch(ensembl_to_entrez(names(ranked_ensembl)), error = function(e) NULL)
  
  if (!is.null(entrez_map) && nrow(entrez_map) >= 100) {
    # Build Entrez-keyed ranked list (keep highest abs rank per Entrez ID)
    rank_df <- data.frame(ENSEMBL = names(ranked_ensembl), rank_score = unname(ranked_ensembl)) |>
      inner_join(entrez_map, by = "ENSEMBL") |>
      arrange(desc(abs(rank_score))) |>
      filter(!duplicated(ENTREZID))
    ranked_entrez <- setNames(rank_df$rank_score, rank_df$ENTREZID)
    ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)
    
    results[["GSEA_KEGG"]] <- tryCatch(
      gseKEGG(geneList = ranked_entrez, organism = "hsa",
              minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
      error = function(e) { message("    GSEA KEGG error: ", e$message); NULL }
    )
    
    # Hallmark (MSigDB): uses Entrez (ncbi_gene)
    hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |>
      select(gs_name, ncbi_gene) |>
      mutate(ncbi_gene = as.character(ncbi_gene))
    
    results[["GSEA_Hallmark"]] <- tryCatch(
      GSEA(geneList = ranked_entrez, TERM2GENE = hallmark,
           minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE),
      error = function(e) { message("    GSEA Hallmark error: ", e$message); NULL }
    )
  } else {
    message("    Entrez conversion yielded too few genes; skipping KEGG/Hallmark GSEA.")
  }
  
  results
}

export_enrichment <- function(result_list, out_dir, prefix) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (name in names(result_list)) {
    res <- result_list[[name]]
    if (is.null(res)) next
    df <- tryCatch(as.data.frame(res), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    export(df, file.path(out_dir, paste0(prefix, "_", name, ".csv")))
  }
}

.pretty_set_label <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^GOBP_", "", x); x <- sub("^GOCC_", "", x); x <- sub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  acronyms <- c("DNA","RNA","TNFA","NFKB","IL1","IL2","IL6","JAK","STAT1","STAT3",
                "STAT5","KRAS","MTORC1","TGF","E2F","UV","WNT","MYC")
  vapply(strsplit(x, " ", fixed = TRUE), function(words) {
    hit <- toupper(words) %in% acronyms
    words[hit] <- toupper(words[hit])
    paste(words, collapse = " ")
  }, character(1))
}

save_dotplot <- function(res, title, path, width = 10, height = 8, show = 20) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  is_gsea <- all(c("NES", "setSize") %in% colnames(df))
  if (is_gsea) {
    df <- df[order(df$p.adjust), , drop = FALSE]
    df <- head(df, show)
    df$Label <- stringr::str_wrap(.pretty_set_label(df$Description), 34)
    p <- ggplot(df, aes(x = NES, y = reorder(Label, NES))) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
      geom_segment(aes(x = 0, xend = NES, yend = reorder(Label, NES)), color = "grey80", linewidth = 0.5) +
      geom_point(aes(size = setSize, fill = p.adjust), shape = 21, color = "grey30") +
      scale_fill_gradient(low = "#b2182b", high = "#2166ac", trans = "log10", name = "p.adjust") +
      scale_size_continuous(range = c(3, 10), name = "Set size") +
      labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
      theme_bw(base_size = 13)
  } else {
    p <- dotplot(res, showCategory = show) + ggtitle(title)
  }
  ggsave(path, plot = p, width = width, height = height, dpi = 600, bg = "white")
}

save_gsea_plot <- function(res, title, path, n = 3) {
  if (is.null(res)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(res), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  ids <- head(df$ID, n)
  p <- gseaplot2(res, geneSetID = ids, title = title)
  ggsave(path, plot = p, width = 12, height = 6 * ceiling(n / 2), dpi = 300)
}

# --- Driver (single condition: ASD vs Healthy Control) -----------------------

message("\n=== ", toupper(PREFIX), " ENRICHMENT ===")

meta_enrich <- import(META_FILE) |> mutate(Gene_ID = as.character(Gene_ID))
meta_enrich <- meta_enrich |> mutate(rem_padj = p.adjust(randomP, method = "BH"))

message("  Genes in meta: ", nrow(meta_enrich))

universe_ensembl <- meta_enrich |> filter(!is.na(Gene_ID)) |> pull(Gene_ID) |> unique()
sig_enrich <- meta_enrich |> filter(rem_padj < PADJ_CUTOFF, abs(randomSummary) >= LFC_CUTOFF)

entrez_up   <- sig_enrich |> filter(randomSummary > 0) |> pull(Gene_ID) |> unique()
entrez_down <- sig_enrich |> filter(randomSummary < 0) |> pull(Gene_ID) |> unique()
entrez_all  <- sig_enrich |> pull(Gene_ID) |> unique()

message("  Significant DEGs (padj<0.05, |LFC|>=1): ", length(entrez_all))
message("  Up: ", length(entrez_up), "  Down: ", length(entrez_down))

message("  Running ORA...")
ora_up   <- run_ora(entrez_up,   universe_ensembl)
ora_down <- run_ora(entrez_down, universe_ensembl)
ora_all  <- run_ora(entrez_all,  universe_ensembl)

export_enrichment(ora_up,   OUT_CSV, paste0(PREFIX, "_ora_up"))
export_enrichment(ora_down, OUT_CSV, paste0(PREFIX, "_ora_down"))
export_enrichment(ora_all,  OUT_CSV, paste0(PREFIX, "_ora_all"))

for (name in names(ora_all)) {
  save_dotplot(ora_all[[name]], paste(PREFIX, name),
               file.path(OUT_FIG, paste0(PREFIX, "_ora_all_", name, "_dotplot.png")))
}

message("  Running GSEA...")
ranked <- build_ranked_list(meta_enrich)
message("  Ranked genes for GSEA: ", length(ranked))
gsea_res <- run_gsea(ranked)
export_enrichment(gsea_res, OUT_CSV, paste0(PREFIX, "_gsea"))

for (name in names(gsea_res)) {
  save_dotplot(gsea_res[[name]], paste(PREFIX, name),
               file.path(OUT_FIG, paste0(PREFIX, "_", name, "_dotplot.png")))
  save_gsea_plot(gsea_res[[name]], paste(PREFIX, name, "- Top Pathways"),
                 file.path(OUT_FIG, paste0(PREFIX, "_", name, "_enrichplot.png")))
}
