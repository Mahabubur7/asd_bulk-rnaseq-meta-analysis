# Meta-Analysis of RNA-seq data using MetaVolcanoR
# Autism Spectrum Disorder iPSC models (ASD vs Healthy Control) across GEO datasets

# Load required packages
library(MetaVolcanoR)
library(metafor)
library(tidyverse)
library(rio)

# Pin dplyr verbs (plyr, if attached, masks these)
mutate <- dplyr::mutate; summarise <- dplyr::summarise; summarize <- dplyr::summarize
arrange <- dplyr::arrange; rename <- dplyr::rename; count <- dplyr::count
desc <- dplyr::desc; select <- dplyr::select; filter <- dplyr::filter


# Output directories
dir.create("results/tables/meta-analysis",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures/meta-analysis", showWarnings = FALSE, recursive = TRUE)

# Compatibility patch: MetaVolcanoR 1.0.1 x ggplot2 4.0
# ggplot2 4.0 plot objects are S7; MetaVolcanoR's S4 "MetaVolcano" class types
# its plot slots as "gg" and rejects them at construction (validObject error).
# Relax the plot slots to "ANY" so the object builds; metaresult is unaffected.
local({
  ns <- asNamespace("MetaVolcanoR")
  nm <- ".__C__MetaVolcano"
  if (bindingIsLocked(nm, ns)) unlockBinding(nm, ns)
  suppressWarnings(setClass("MetaVolcano",
                            representation(input = "data.frame", inputnames = "character",
                                           metaresult = "data.frame", MetaVolcano = "ANY", degfreq = "ANY"),
                            where = ns))
})

# Gene annotation lookup, reused from the pre-annotated tables produced by
# 02_annotation.R (results/tables/annotated). No re-annotation / network calls.
anno_files <- list.files("results/tables/annotated", pattern = "[.]csv$", full.names = TRUE)
gene_annotation <- anno_files |>
  lapply(function(f) {
    import(f) |> select(any_of(c("Gene_ID", "Gene_Symbol", "Gene_Description")))
  }) |>
  bind_rows() |>
  mutate(Gene_ID = as.character(Gene_ID)) |>
  distinct(Gene_ID, .keep_all = TRUE)

annotate_genes <- function(df, id_col = "Gene_ID") {
  df |>
    mutate(!!id_col := as.character(.data[[id_col]])) |>
    left_join(gene_annotation, by = id_col)
}

# Data pre-processing
# Auto-detect datasets from DESeq2 output
deg_files <- list.files("results/tables/DESeq2", pattern = "\\.csv$", full.names = TRUE)
geo_ids   <- tools::file_path_sans_ext(basename(deg_files))

# Read all studies into a named list
studies <- lapply(deg_files, import)
names(studies) <- geo_ids

# --- Manual fallback implementations -----------------------------------
# The iza-mcac/MetaVolcanoR fork's rem_mv() has a confirmed bug on this
# dataset. These functions reproduce the same statistics directly via
# metafor::rma() (REM) and Fisher's method (Combining), and are used
# automatically if the official functions fail.

run_manual_rem <- function(studies, genenamecol = "Gene_ID",
                            foldchangecol = "log2FoldChange", vcol = "lfcSE") {
  prepped <- Map(function(df, i) {
    df <- df |> select(all_of(c(genenamecol, foldchangecol, vcol)))
    df$vi <- df[[vcol]]^2
    df <- df |> select(all_of(genenamecol), all_of(foldchangecol), vi)
    names(df)[names(df) == foldchangecol] <- paste0("fc_", i)
    names(df)[names(df) == "vi"] <- paste0("vi_", i)
    df
  }, studies, seq_along(studies))

  merged <- Reduce(function(x, y) merge(x, y, by = genenamecol, all = TRUE), prepped)
  fc_cols <- grep("^fc_", names(merged), value = TRUE)
  vi_cols <- grep("^vi_", names(merged), value = TRUE)

  fit_one <- function(row) {
    fc <- as.numeric(row[fc_cols]); vi <- as.numeric(row[vi_cols])
    keep <- !is.na(fc) & !is.na(vi)
    fc <- fc[keep]; vi <- vi[keep]
    if (length(fc) < 2) {
      return(c(randomSummary = NA_real_, randomP = NA_real_,
               randomCi.lb = NA_real_, randomCi.ub = NA_real_))
    }
    fit <- tryCatch(metafor::rma(yi = fc, vi = vi, method = "REML"),
                     error = function(e) NULL)
    if (is.null(fit)) {
      return(c(randomSummary = NA_real_, randomP = NA_real_,
               randomCi.lb = NA_real_, randomCi.ub = NA_real_))
    }
    c(randomSummary = as.numeric(fit$beta), randomP = fit$pval,
      randomCi.lb = fit$ci.lb, randomCi.ub = fit$ci.ub)
  }

  rem_stats <- t(apply(merged, 1, fit_one))
  cbind(merged[genenamecol], as.data.frame(rem_stats))
}

run_manual_combining <- function(studies, genenamecol = "Gene_ID",
                                  foldchangecol = "log2FoldChange", pcol = "padj") {
  prepped <- Map(function(df, i) {
    df <- df |> select(all_of(c(genenamecol, foldchangecol, pcol)))
    names(df)[names(df) == foldchangecol] <- paste0("fc_", i)
    names(df)[names(df) == pcol] <- paste0("p_", i)
    df
  }, studies, seq_along(studies))

  merged <- Reduce(function(x, y) merge(x, y, by = genenamecol, all = TRUE), prepped)
  fc_cols <- grep("^fc_", names(merged), value = TRUE)
  p_cols  <- grep("^p_",  names(merged), value = TRUE)

  fisher_combine <- function(p) {
    p <- p[!is.na(p)]
    if (length(p) < 2) return(NA_real_)
    p[p <= 0] <- .Machine$double.eps
    stat <- -2 * sum(log(p))
    pchisq(stat, df = 2 * length(p), lower.tail = FALSE)
  }

  merged$metafc <- rowMeans(merged[fc_cols], na.rm = TRUE)
  merged$metap  <- apply(merged[p_cols], 1, fisher_combine)
  merged |> select(all_of(genenamecol), metafc, metap)
}

# 04. Meta-Analysis -- Random Effect Model
meta_results <- tryCatch({
  meta_degs_rem <- rem_mv(
    diffexp       = studies,
    pcriteria     = "padj",
    foldchangecol = "log2FoldChange",
    genenamecol   = "Gene_ID",
    geneidcol     = NULL,
    collaps       = TRUE,
    vcol          = "lfcSE",
    cvar          = FALSE,
    metathr       = 0.01,
    jobname       = "MetaVolcano_REM_ASD",
    outputfolder  = "results/figures/meta-analysis/",
    draw          = "PDF",
    ncores        = 1 # adjust as per your PC core! 
  )
  meta_degs_rem@metaresult
}, error = function(e) {
  message("rem_mv() failed (", conditionMessage(e), ") -- falling back to a manual ",
          "random-effects meta-analysis via metafor::rma(), which reproduces the ",
          "same statistics and is confirmed to work on this dataset.")
  run_manual_rem(studies, genenamecol = "Gene_ID",
                  foldchangecol = "log2FoldChange", vcol = "lfcSE")
})

# Annotate + export full REM result
annotated_results <- annotate_genes(meta_results)
export(annotated_results, "results/tables/meta-analysis/random_effect_model.csv")

# Filter DEGs: significant (randomP < 0.05) and reliable effect (|randomSummary| >= 1)
key_genes <- annotated_results |>
  filter(randomP < 0.05, abs(randomSummary) >= 1)
export(key_genes, "results/tables/meta-analysis/filtered_meta_degs.csv")

# Meta-Analysis -- Combining approach (Mean)
combined_meta <- tryCatch({
  meta_degs_comb <- combining_mv(
    diffexp       = studies,
    pcriteria     = "padj",
    foldchangecol = "log2FoldChange",
    genenamecol   = "Gene_ID",
    metafc        = "Mean",
    metathr       = 0.01,
    collaps       = TRUE,
    jobname       = "MetaVolcano_Combining_ASD",
    outputfolder  = "results/figures/meta-analysis/",
    draw          = "PDF"
  )
  meta_degs_comb@metaresult
}, error = function(e) {
  message("combining_mv() failed (", conditionMessage(e), ") -- falling back to a ",
          "manual combining approach (mean log2FC + Fisher's method for combined ",
          "p-values across datasets).")
  run_manual_combining(studies, genenamecol = "Gene_ID",
                        foldchangecol = "log2FoldChange", pcol = "padj")
})
export(annotate_genes(combined_meta),
       "results/tables/meta-analysis/meta_combining_mean.csv")

# Cross-reference with Venn intersection (genes DE across all datasets)
venn_file <- "results/tables/venn/common_genes_all_datasets.csv"
if (file.exists(venn_file)) {
  intersect_genes <- import(venn_file)
  
  intersects <- combined_meta |>
    filter(Gene_ID %in% as.character(intersect_genes$Gene_ID))
  
  intersec_results <- annotate_genes(intersects)
  export(intersec_results,
         "results/tables/meta-analysis/meta_combining_mean_intersect_genes.csv")
  message("Intersection genes annotated: ", nrow(intersec_results))
} else {
  message("Venn intersection file not found; skipping intersect cross-reference. Run 03_venn.R first.")
}
