#!/usr/bin/env Rscript
## =============================================================================
## 01_load_cellranger.R — CellRanger matrices -> SingleCellExperiment
##
## Generic version: Samples are either auto-discovered from the project directory or
## supplied explicitly, and the directory layout is described by a template rather than baked into the loop.
## 
##
## CONFIGURE THREE THINGS
##   PROJECT   Data location
##   SAMPLES   NULL to auto-discover, or an explicit character vector
##   LAYOUT    where the matrix sits under each sample directory
##
## LAYOUT is a template with {sample} substituted per sample. Common cases:
##   "{sample}/cellranger_count/{sample}/outs/filtered_feature_bc_matrix"  <- nested
##   "{sample}/outs/filtered_feature_bc_matrix"                            <- standard
##   "{sample}/filtered_feature_bc_matrix"                                 <- flat
##   "{sample}"                                                            <- GEO-style
## 
## =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(DropletUtils)
})

## ---- Config -----------------------------------------------------------------
PROJECT <- "/g/data/pq84/single_cell/10X/scRNA_auto/all_data/GIMR_GGP_CIRCA_20250513_5P"

## NULL = auto-discover every directory under PROJECT that contains a matrix.
## Or set explicitly, e.g. c("Pool_A","Pool_B") or paste0("Pool_", LETTERS[1:4]).
SAMPLES <- NULL

## Candidate layouts, tried in order. First hit wins, per sample.
LAYOUTS <- c(
  "{sample}/cellranger_count/{sample}/outs/filtered_feature_bc_matrix",
  "{sample}/outs/filtered_feature_bc_matrix",
  "{sample}/filtered_feature_bc_matrix",
  "{sample}"
)

MATRIX_TYPE <- "filtered"   # "filtered" or "raw"; swaps the layout templates
OVERWRITE   <- FALSE        # TRUE to re-load samples that already have an .rds
COMPRESS    <- "xz"         # "xz" (small, slow) / "gzip" (fast) / FALSE

OUT_SCE <- file.path(PROJECT, "data", "SCEs")
OUT_BC  <- file.path(PROJECT, "data", "barcodes")
OUT_FIG <- file.path(PROJECT, "output", "figures", "QC")
for (d in c(OUT_SCE, OUT_BC, OUT_FIG))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (MATRIX_TYPE == "raw")
  LAYOUTS <- gsub("filtered_feature_bc_matrix", "raw_feature_bc_matrix", LAYOUTS)

## ---- Path resolution ---------------------------------------------------------
is_matrix_dir <- function(p) {
  dir.exists(p) &&
    any(grepl("^(matrix\\.mtx|barcodes\\.tsv|features\\.tsv|genes\\.tsv)",
              list.files(p)))
}

resolve_sample <- function(sample, project = PROJECT, layouts = LAYOUTS) {
  for (tpl in layouts) {
    ## plain substitution — avoids glue's .envir requirements entirely
    p <- file.path(project, gsub("{sample}", sample, tpl, fixed = TRUE))
    if (is_matrix_dir(p)) return(p)
  }
  NA_character_
}

## ---- Sample discovery --------------------------------------------------------
if (is.null(SAMPLES)) {
  cand <- list.dirs(PROJECT, recursive = FALSE, full.names = FALSE)
  cand <- setdiff(cand, c("data", "output", "code", "scripts", "logs", "figures"))
  SAMPLES <- cand[!is.na(vapply(cand, resolve_sample, character(1)))]
  if (!length(SAMPLES))
    stop("No samples auto-discovered under '", PROJECT, "'.\n",
         "  Set SAMPLES explicitly, or add the correct path to LAYOUTS.\n",
         "  Directories seen: ", paste(head(cand, 20), collapse = ", "))
  message("auto-discovered ", length(SAMPLES), " sample(s)")
}

## ---- Pre-flight: show every resolved path before loading anything ------------
paths <- vapply(SAMPLES, resolve_sample, character(1))
cat("\n=== Path check ===\n")
for (s in SAMPLES)
  cat(sprintf("  %-28s %s\n    %s\n", s,
              if (is.na(paths[[s]])) "NOT FOUND" else "ok",
              if (is.na(paths[[s]])) "-" else sub(PROJECT, ".", paths[[s]], fixed = TRUE)))
if (anyNA(paths))
  stop(sum(is.na(paths)), " sample(s) unresolved — fix LAYOUTS or SAMPLES ",
       "before running.")

## ---- Load --------------------------------------------------------------------
summary_rows <- list()

for (s in SAMPLES) {
  out_file <- file.path(OUT_SCE, paste0(s, ".CellRanger.SCE.rds"))

  if (file.exists(out_file) && !OVERWRITE) {
    cat("\n", s, ": already loaded, skipping\n", sep = "")
    next
  }

  cat("\n=== ", s, " ===\n", sep = "")
  sce <- read10xCounts(samples = setNames(paths[[s]], s), col.names = TRUE)
  stopifnot(!anyDuplicated(colnames(sce)))

  ## carry the sample id in colData so downstream merges can group on it
  sce$sample_id <- s

  ## split antibody capture / CRISPR features into altExps if present
  if ("Type" %in% colnames(rowData(sce)) &&
      length(unique(rowData(sce)$Type)) > 1) {
    sce <- splitAltExps(sce, rowData(sce)$Type, "Gene Expression")
    cat("  altExps: ", paste(altExpNames(sce), collapse = ", "), "\n", sep = "")
  }

  cat(sprintf("  %d cells x %d genes\n", ncol(sce), nrow(sce)))
  saveRDS(sce, out_file, compress = COMPRESS)
  writeLines(colnames(sce), file.path(OUT_BC, paste0(s, "_barcodes.txt")))
  cat("  saved: ", basename(out_file), "\n", sep = "")

  summary_rows[[s]] <- data.frame(
    sample = s, cells = ncol(sce), genes = nrow(sce),
    matrix = MATRIX_TYPE, stringsAsFactors = FALSE)
}

## ---- Summary -----------------------------------------------------------------
if (length(summary_rows)) {
  sm <- do.call(rbind, summary_rows)
  cat("\n=== LOAD SUMMARY ===\n"); print(sm, row.names = FALSE)
  cat("total cells:", sum(sm$cells), "\n")
  utils::write.csv(sm, file.path(PROJECT, "output", "load_summary.csv"),
                   row.names = FALSE)
} else {
  cat("\nNothing loaded (all samples already present; set OVERWRITE = TRUE to redo)\n")
}
cat("SCEs in ", OUT_SCE, "\n", sep = "")
