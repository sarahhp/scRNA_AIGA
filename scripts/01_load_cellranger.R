# ══════════════════════════════════════════════════════════════════════════════
# 01_load_cellranger.R — Load CellRanger filtered matrices → SCE
# ══════════════════════════════════════════════════════════════════════════════

setwd("/g/data/pq84/single_cell/10X/scRNA_auto/all_data/GIMR_GGP_CIRCA_20250513_5P/")

library(here)
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(glue)

here()

# Verify paths exist before the loop
for (i in c("A", "B")) {
  p <- here(paste0("Pool_", i), "cellranger_count",
            paste0("Pool_", i), "outs", "filtered_feature_bc_matrix")
  cat(i, ":", ifelse(dir.exists(p), "✔ EXISTS", "✘ NOT FOUND"), "\n")
  cat("  ", p, "\n")
}

dir.create(here("data", "SCEs"),     recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "barcodes"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures", "QC"), recursive = TRUE, showWarnings = FALSE)

for (i in c("A", "B")) {

  out_file <- here("data", "SCEs", glue("Pool_{i}.CellRanger.SCE.rds"))

  if (file.exists(out_file)) {
    cat("✔ Pool", i, "already loaded — skipping\n")
    next
  }

  cat("\n── Loading Pool", i, "──────────────────────────\n")

  capture_name <- paste0("Pool_", i)

  # Path to CellRanger filtered matrix
  filtered_path <- here(
    capture_name, "cellranger_count", capture_name,
    "outs", "filtered_feature_bc_matrix"
  )

  sce <- read10xCounts(
    samples   = setNames(filtered_path, capture_name),
    col.names = TRUE
  )
  stopifnot(!anyDuplicated(colnames(sce)))

  # Split alternative experiments if present
  if ("Type" %in% colnames(rowData(sce))) {
    sce <- splitAltExps(sce, rowData(sce)$Type, "Gene Expression")
  }

  cat("  Cells loaded:", ncol(sce), "\n")
  cat("  Genes loaded:", nrow(sce), "\n")

  saveRDS(sce, out_file, compress = "xz")
  cat("  ✔ Saved:", out_file, "\n")

  # Save barcodes
  writeLines(
    colnames(sce),
    here("data", "barcodes", glue("Pool_{i}_barcodes.txt"))
  )
}
