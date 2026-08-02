# ══════════════════════════════════════════════════════════════════════════════
# 03_preprocess.R
# QC, gene annotation, doublet removal, vireo assignment
# Pools A and B
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Fix here() root before anything else ────────────────────────────────────
PROJECT_ROOT <- "/g/data/pq84/single_cell/10X/scRNA_auto/all_data/GIMR_GGP_CIRCA_20250513_5P"

setwd(PROJECT_ROOT)
try(detach("package:here", unload = TRUE), silent = TRUE)
library(here)

if (here() != PROJECT_ROOT) {
  stop(glue::glue(
    "here() still wrong: {here()}\n",
    "Start R from: {PROJECT_ROOT}"
  ))
}
cat("✔ Project root:", here(), "\n")

# ── 1. Libraries ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(scater)
  library(scuttle)
  library(EnsDb.Hsapiens.v86)
  library(msigdbr)
  library(glue)
  library(dplyr)
  library(DropletUtils)
  library(scDblFinder)
})

dir.create(here("output", "figures", "QC"), recursive = TRUE, showWarnings = FALSE)

# ── 2. Load central donor metadata ─────────────────────────────────────────────
donor_meta_file <- here("data", "meta_data", "donor_meta.csv")

if (!file.exists(donor_meta_file)) {
  stop(glue(
    "donor_meta.csv not found: {donor_meta_file}\n",
    "Run 00_create_donor_metadata.R first."
  ))
}

donor_meta <- read.csv(donor_meta_file, stringsAsFactors = FALSE)

cat("✔ Donor metadata loaded:\n")
print(donor_meta[, c("donor_short", "cohort", "disease", "pool")])

# Build vireo long ID → short name map from CSV (no hard-coded IDs)
donor_mapping <- setNames(
  c(donor_meta$donor_short, "doublet", "unassigned"),
  c(donor_meta$vireo_id,    "doublet", "unassigned")
)

donor_levels <- c(donor_meta$donor_short, "doublet", "unassigned")

# ── 3. Ribosomal gene set (computed once) ──────────────────────────────────────
c2_sets  <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_set <- character(0)

# ══════════════════════════════════════════════════════════════════════════════
# PREPROCESS EACH POOL
# ══════════════════════════════════════════════════════════════════════════════
for (i in c("A", "B")) {
  
  out_file <- here("data", "SCEs", glue("Pool_{i}.preprocessed.SCE.rds"))
  
  if (file.exists(out_file)) {
    cat(glue("✔ Pool {i} already preprocessed — skipping\n"))
    next
  }
  
  cat(glue("\n{strrep('─', 50)}\n  Preprocessing Pool {i}\n{strrep('─', 50)}\n"))
  
  # ── Load CellRanger SCE ─────────────────────────────────────────────────────
  sce_file <- here("data", "SCEs", glue("Pool_{i}.CellRanger.SCE.rds"))
  if (!file.exists(sce_file)) stop(glue("CellRanger SCE not found: {sce_file}"))
  
  sce <- readRDS(sce_file)
  cat("  Cells loaded:", ncol(sce), "\n")
  
  # ── Identify empty droplets ─────────────────────────────────────────────────
  sce$capture <- factor(sce$Sample)
  capture_names <- setNames(levels(sce$capture), levels(sce$capture))
  
  empties <- do.call(rbind, lapply(capture_names, function(cn) {
    emp_file <- here("data", "emptyDrops", paste0(cn, ".emptyDrops.rds"))
    if (!file.exists(emp_file)) {
      cat(glue("  ⚠ emptyDrops file not found for {cn} — skipping empty drop filter\n"))
      return(NULL)
    }
    empties <- readRDS(emp_file)
    empties$capture <- cn
    empties
  }))
  
  if (!is.null(empties)) {
    cat("  Non-empty droplets per capture:\n")
    print(tapply(empties$FDR, empties$capture,
                 function(x) sum(x <= 0.001, na.rm = TRUE)))
    # Remove empty droplets
    sce <- sce[, which(empties$FDR <= 0.001)]
    cat("  Cells after empty drop removal:", ncol(sce), "\n")
  }
  
  # ── Doublet detection ───────────────────────────────────────────────────────
  dbl_out <- here("data", "SCEs", glue("Pool_{i}_doublets.rds"))
  
  if (!file.exists(dbl_out)) {
    cat("  Running doublet detection...\n")
    sce <- scDblFinder(sce)
    dbl_df <- data.frame(colData(sce)) %>%
      tibble::rownames_to_column(var = "cell")
    saveRDS(dbl_df, file = dbl_out)
    cat("  ✔ Doublet calls saved\n")
  } else {
    cat("  ✔ Doublet calls already exist\n")
    sce <- scDblFinder(sce)
  }
  
  cat("  Doublet summary:\n")
  print(table(sce$scDblFinder.class))
  
  # ── Gene annotation ─────────────────────────────────────────────────────────
  rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$ID,
    rowData(sce)$Symbol
  )
  
  location <- mapIds(
    x       = EnsDb.Hsapiens.v86,
    keys    = rowData(sce)$ID,
    keytype = "GENEID",
    column  = "SEQNAME"
  )
  rowData(sce)$CHR <- location
  
  biotype <- mapIds(
    x       = EnsDb.Hsapiens.v86,
    keys    = rowData(sce)$ID,
    keytype = "GENEID",
    column  = "GENEBIOTYPE"
  )
  rowData(sce)$ENSEMBL.GENEBIOTYPE <- biotype
  
  # ── Gene sets ───────────────────────────────────────────────────────────────
  mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
  
  ribo_genes <- grep("^RP(S|L)", rownames(sce), value = TRUE)
  ribo_set   <- union(
    ribo_genes,
    c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_gene_symbol
  )
  
  sex_set        <- rownames(sce)[rowData(sce)$CHR %in% c("X", "Y")]
  pseudogene_set <- rownames(sce)[grepl("pseudogene",
                                        rowData(sce)$ENSEMBL.GENEBIOTYPE,
                                        ignore.case = TRUE)]
  
  cat(glue("  Mito genes: {length(mito_set)} | Ribo genes: {length(ribo_set)}\n"))
  cat(glue("  Sex genes: {length(sex_set)} | Pseudogenes: {length(pseudogene_set)}\n"))
  
  # ── Per-cell QC metrics ─────────────────────────────────────────────────────
  sce <- addPerCellQC(
    sce,
    subsets = list(
      Mito = which(rownames(sce) %in% mito_set),
      Ribo = which(rownames(sce) %in% ribo_set)
    )
  )
  
  # ── QC plots ────────────────────────────────────────────────────────────────
  pdf(here("output", "figures", "QC", glue("Pool_{i}_QC.pdf")),
      width = 12, height = 8)
  print(
    plotColData(sce, "sum",                  colour_by = "sum")      +
      scale_y_log10() +
      plotColData(sce, "detected",             colour_by = "detected") +
      plotColData(sce, "subsets_Mito_percent", colour_by = "subsets_Mito_percent") +
      plotColData(sce, "subsets_Ribo_percent", colour_by = "subsets_Ribo_percent")
  )
  dev.off()
  cat("  ✔ QC plots saved\n")
  
  # ── Filter high mitochondrial cells ─────────────────────────────────────────
  mito_drop <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher")
  cat("  Cells removed (high mito):", sum(mito_drop), "\n")
  sce <- sce[, !mito_drop]
  
  # ── Remove doublets ─────────────────────────────────────────────────────────
  sce <- sce[, sce$scDblFinder.class == "singlet"]
  cat("  Cells after doublet removal:", ncol(sce), "\n")
  
  # ── Assign vireo donor IDs ───────────────────────────────────────────────────
  vireo_file <- here(glue("vireo_Flt_Pool_{i}_out"), "donor_ids.tsv")
  if (!file.exists(vireo_file)) stop(glue("Vireo file not found: {vireo_file}"))
  
  vireo_df <- read.table(vireo_file, header = TRUE)
  
  j <- match(colnames(sce), vireo_df$cell)
  
  # Use donor_mapping from donor_meta.csv — no hard-coded vireo IDs
  sce$genetic_donor <- factor(
    dplyr::recode(vireo_df$donor_id[j], !!!donor_mapping),
    levels = donor_levels
  )
  
  # ── Attach donor metadata from CSV ──────────────────────────────────────────
  meta_join <- donor_meta[
    match(as.character(sce$genetic_donor), donor_meta$donor_short),
    c("cohort", "disease", "pool", "sex")
  ]
  sce$cohort  <- meta_join$cohort
  sce$disease <- meta_join$disease
  sce$sex     <- meta_join$sex
  
  # ── Mark cell status ─────────────────────────────────────────────────────────
  sce$cell_status <- "cell"
  
  cat("  Donor assignment:\n")
  print(table(sce$genetic_donor))
  cat("  Cohort assignment:\n")
  print(table(sce$cohort, useNA = "ifany"))
  
  # ── Save ─────────────────────────────────────────────────────────────────────
  saveRDS(sce, out_file, compress = "xz")
  cat(glue("  ✔ Saved: {out_file}\n"))
}

cat(glue("\n{strrep('═', 50)}\n  PREPROCESSING COMPLETE\n{strrep('═', 50)}\n"))

