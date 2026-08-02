# ══════════════════════════════════════════════════════════════════════════════
# 04_integrate_annotate.R
# Filter, combine, integrate (SCTransform), annotate with Azimuth
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
  library(scuttle)
  library(Seurat)
  library(Azimuth)
  library(msigdbr)
  library(dplyr)
  library(glue)
  library(patchwork)
  library(ggplot2)
  library(scales)
})

dir.create(here("output", "figures", "Integration"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("output", "figures", "CellType"),    recursive = TRUE, showWarnings = FALSE)

# ── 2. Load central donor metadata ─────────────────────────────────────────────
donor_meta_file <- here("data", "metadata", "donor_meta.csv")

if (!file.exists(donor_meta_file)) {
  stop(glue(
    "donor_meta.csv not found: {donor_meta_file}\n",
    "Run 00_create_donor_metadata.R first."
  ))
}

donor_meta <- read.csv(donor_meta_file, stringsAsFactors = FALSE)

# All derived from CSV — nothing hard-coded
cohort_levels  <- c("Healthy",
                    unique(donor_meta$cohort[donor_meta$cohort != "Healthy"]))
healthy_donors <- donor_meta$donor_short[donor_meta$cohort == "Healthy"]
donor_levels   <- c(donor_meta$donor_short, "doublet", "unassigned")

cat("✔ Donor metadata loaded:\n")
print(donor_meta[, c("donor_short", "cohort", "disease", "pool")])
cat(glue("Cohort levels: {paste(cohort_levels, collapse = ' > ')}\n"))

# ── Helper: attach cohort / disease / sex to Seurat object ────────────────────
attach_donor_meta <- function(seu_obj, meta = donor_meta) {
  idx             <- match(as.character(seu_obj$genetic_donor), meta$donor_short)
  seu_obj$cohort  <- factor(meta$cohort[idx], levels = cohort_levels)
  seu_obj$disease <- meta$disease[idx]
  seu_obj$sex     <- meta$sex[idx]
  seu_obj
}

# ── 3. Gene filter function ────────────────────────────────────────────────────
c2_sets <- msigdbr(species = "Homo sapiens", category = "C2")

filter_genes <- function(sce) {
  mito   <- rownames(sce)[rowData(sce)$CHR == "MT"]
  ribo   <- union(
    grep("^RP(S|L)", rownames(sce), value = TRUE),
    c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", ]$human_gene_symbol
  )
  sex    <- rownames(sce)[rowData(sce)$CHR %in% c("X", "Y")]
  pseudo <- rownames(sce)[grepl("pseudogene",
                                rowData(sce)$ENSEMBL.GENEBIOTYPE,
                                ignore.case = TRUE)]
  
  uninformative <- rownames(sce) %in% c(mito, ribo, sex, pseudo)
  junk          <- sce$genetic_donor %in% c("doublet", "unassigned") |
    sce$cell_status   == "empty_droplet"
  
  cat(glue("  Uninformative genes removed: {sum(uninformative)}\n"))
  cat(glue("  Junk cells removed:          {sum(junk)}\n"))
  
  sce[!uninformative, !junk]
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN — build or load annotated Seurat object
# ══════════════════════════════════════════════════════════════════════════════
out_file <- here("data", "SCEs", "Pools_A_B.annotated.SEU.rds")

if (file.exists(out_file)) {
  cat("✔ Annotated object exists — loading\n")
  seu <- readRDS(out_file)
  
} else {
  
  # ── Load and filter Pool A ───────────────────────────────────────────────────
  cat(glue("\n{strrep('─', 50)}\n  Filtering Pool A\n{strrep('─', 50)}\n"))
  sce_A_file <- here("data", "SCEs", "Pool_A.preprocessed.SCE.rds")
  if (!file.exists(sce_A_file)) stop(glue("Missing: {sce_A_file}"))
  
  sce_A    <- readRDS(sce_A_file)
  sceFlt_A <- filter_genes(sce_A)
  sceFlt_A$pool <- "Pool_A"
  colnames(sceFlt_A) <- paste0(colnames(sceFlt_A), "_A")
  cat(glue("  Pool A cells after filter: {ncol(sceFlt_A)}\n"))
  
  # ── Load and filter Pool B ───────────────────────────────────────────────────
  cat(glue("\n{strrep('─', 50)}\n  Filtering Pool B\n{strrep('─', 50)}\n"))
  sce_B_file <- here("data", "SCEs", "Pool_B.preprocessed.SCE.rds")
  if (!file.exists(sce_B_file)) stop(glue("Missing: {sce_B_file}"))
  
  sce_B    <- readRDS(sce_B_file)
  sceFlt_B <- filter_genes(sce_B)
  sceFlt_B$pool <- "Pool_B"
  colnames(sceFlt_B) <- paste0(colnames(sceFlt_B), "_B")
  cat(glue("  Pool B cells after filter: {ncol(sceFlt_B)}\n"))
  
  # ── Shared genes ─────────────────────────────────────────────────────────────
  shared_genes <- intersect(rownames(sceFlt_A), rownames(sceFlt_B))
  cat(glue("\n  Shared genes: {length(shared_genes)}\n"))
  sceFlt_A <- sceFlt_A[shared_genes, ]
  sceFlt_B <- sceFlt_B[shared_genes, ]
  
  # ── Remove low-abundance genes ───────────────────────────────────────────────
  combined_counts <- cbind(counts(sceFlt_A), counts(sceFlt_B))
  keep            <- rowSums(combined_counts > 0) > 20
  cat(glue("  Genes kept (>20 cells expressing): {sum(keep)}\n"))
  sceFlt_A <- sceFlt_A[keep, ]
  sceFlt_B <- sceFlt_B[keep, ]
  

  # ── Combine into single SCE ──────────────────────────────────────────────────
  combined_counts  <- cbind(counts(sceFlt_A),  counts(sceFlt_B))
  combined_colData <- rbind(colData(sceFlt_A), colData(sceFlt_B))
  
  sce_combined <- SingleCellExperiment(
    list(counts = combined_counts),
    colData = combined_colData,
    rowData = rowData(sceFlt_A)
  )
  cat(glue("\n  Combined: {ncol(sce_combined)} cells × {nrow(sce_combined)} genes\n"))
  
  # ── Convert to Seurat ────────────────────────────────────────────────────────
  cnt           <- counts(sce_combined)
  rownames(cnt) <- make.unique(rowData(sce_combined)$Symbol)
  
  seu_combined <- CreateSeuratObject(
    counts    = cnt,
    meta.data = data.frame(colData(sce_combined))
  )
  
  # ── SCTransform integration across donors ─────────────────────────────────────
  cat(glue("\n{strrep('─', 50)}\n  Integrating across donors\n{strrep('─', 50)}\n"))
  
  utility_file <- here("code", "utility.R")
  if (!file.exists(utility_file)) stop(glue("utility.R not found: {utility_file}"))
  source(utility_file)
  
  seu_integrated <- intDat(seu_combined, split = "genetic_donor", type = "RNA")
  
  # ── Dimensionality reduction ─────────────────────────────────────────────────
  DefaultAssay(seu_integrated) <- "integrated"
  seu_integrated <- seu_integrated %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)
  
  VariableFeatures(seu_integrated) <- rownames(seu_integrated@assays$integrated)
  
  cat(glue("  PCA and UMAP complete\n"))
  
  # ── Azimuth annotation ───────────────────────────────────────────────────────
  cat(glue("\n{strrep('─', 50)}\n  Running Azimuth annotation\n{strrep('─', 50)}\n"))
  seu <- RunAzimuth(seu_integrated, reference = "pbmcref")
  
  # ── Attach cohort / disease / sex from central CSV ────────────────────────────
  seu <- attach_donor_meta(seu)
  
  # Verify metadata attached
  if (all(is.na(seu$cohort))) {
    stop("cohort is all NA — check donor_short in donor_meta.csv matches seu$genetic_donor")
  }
  
  cat("\n  Cohort summary:\n")
  print(table(seu$cohort))
  cat("\n  Cell type summary (l2):\n")
  print(table(seu$predicted.celltype.l2))
  
  # ── Save ─────────────────────────────────────────────────────────────────────
  saveRDS(seu, out_file)
  cat(glue("  ✔ Saved: {out_file}\n"))
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════
cat(glue("\n{strrep('─', 50)}\n  Generating plots\n{strrep('─', 50)}\n"))

# Cohort colours derived from levels — no hard-coding
cohort_colours <- setNames(
  c("#4DAF4A", "#E41A1C", "#FF7F00", "#984EA3", "#377EB8")[seq_along(cohort_levels)],
  cohort_levels
)

# ── UMAP: integration QC ──────────────────────────────────────────────────────
pdf(here("output", "figures", "Integration", "UMAP_integration.pdf"),
    width = 18, height = 6)
print(
  DimPlot(seu, group.by = "pool",          pt.size = 0.3)                 +
    DimPlot(seu, group.by = "genetic_donor", pt.size = 0.3)                 +
    DimPlot(seu, group.by = "cohort",        pt.size = 0.3,
            cols = cohort_colours)
)
dev.off()
cat("  ✔ UMAP_integration.pdf\n")

# ── UMAP: cell types ──────────────────────────────────────────────────────────
pdf(here("output", "figures", "CellType", "UMAP_celltypes_l1.pdf"),
    width = 10, height = 8)
print(
  DimPlot(seu, group.by = "predicted.celltype.l1",
          label = TRUE, pt.size = 0.3, repel = TRUE) + NoLegend()
)
dev.off()

pdf(here("output", "figures", "CellType", "UMAP_celltypes_l2.pdf"),
    width = 12, height = 10)
print(
  DimPlot(seu, group.by = "predicted.celltype.l2",
          label = TRUE, pt.size = 0.3, repel = TRUE) + NoLegend()
)
dev.off()
cat("  ✔ UMAP_celltypes.pdf\n")

# ── UMAP: per donor ───────────────────────────────────────────────────────────
pdf(here("output", "figures", "CellType", "UMAP_per_donor.pdf"),
    width = 18, height = 10)
print(
  DimPlot(seu, group.by = "predicted.celltype.l2",
          split.by = "genetic_donor",
          label = TRUE, pt.size = 0.3, repel = TRUE, ncol = 3) +
    NoLegend()
)
dev.off()
cat("  ✔ UMAP_per_donor.pdf\n")

# ── Cell type × cohort composition barplot ─────────────────────────────────────
celltype_cohort <- as.data.frame(table(
  CellType = seu$predicted.celltype.l2,
  Cohort   = seu$cohort
)) %>%
  group_by(Cohort) %>%
  mutate(pct = 100 * Freq / sum(Freq)) %>%
  ungroup() %>%
  filter(Freq > 0)

p_comp <- ggplot(celltype_cohort,
                 aes(x = Cohort, y = pct, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  scale_x_discrete(limits = cohort_levels) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(x = NULL, y = "% of cells", fill = "Cell type (l2)",
       title = "Cell type composition per cohort") +
  theme_minimal(base_size = 12) +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        axis.text.x     = element_text(angle = 30, hjust = 1),
        legend.position = "right",
        legend.text     = element_text(size = 7))

ggsave(here("output", "figures", "CellType", "composition_by_cohort.pdf"),
       p_comp, width = 10, height = 7)
cat("  ✔ composition_by_cohort.pdf\n")

# ── Summary tables ────────────────────────────────────────────────────────────
cat("\n  Final cell counts:\n")
cat(glue("  Total cells: {ncol(seu)}\n"))
cat("\n  Cell type (l2) × cohort:\n")
print(table(seu$predicted.celltype.l2, seu$cohort))
cat("\n  Donor × cohort:\n")
print(table(Donor = seu$genetic_donor, Cohort = seu$cohort))

cat(glue("\n{strrep('═', 50)}\n  INTEGRATION & ANNOTATION COMPLETE\n{strrep('═', 50)}\n"))
cat(glue("  Output: {out_file}\n"))

