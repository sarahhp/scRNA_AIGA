# ══════════════════════════════════════════════════════════════════════════════
# 03b_QC_plots.R
# Comprehensive QC visualisation for Pool A and Pool B
# Run AFTER 03_preprocess.R
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Fix here() root ─────────────────────────────────────────────────────────
PROJECT_ROOT <- "/g/data/pq84/single_cell/10X/scRNA_auto/all_data/GIMR_GGP_CIRCA_20250513_5P"

setwd(PROJECT_ROOT)
try(detach("package:here", unload = TRUE), silent = TRUE)
library(here)

if (here() != PROJECT_ROOT) stop(paste("Wrong root:", here()))
cat("✔ Project root:", here(), "\n")

# ── 1. Libraries ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(scater)
  library(scuttle)
  library(DropletUtils)
  library(EnsDb.Hsapiens.v86)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(glue)
})

dir.create(here("output", "figures", "QC"), recursive = TRUE, showWarnings = FALSE)

# ── 2. Load donor metadata ─────────────────────────────────────────────────────
donor_meta_file <- here("data", "meta_data", "donor_meta.csv")
if (!file.exists(donor_meta_file)) {
  stop(glue("donor_meta.csv not found: {donor_meta_file}\nRun 00_create_donor_metadata.R first."))
}
donor_meta <- read.csv(donor_meta_file, stringsAsFactors = FALSE)

cat("✔ Donor metadata loaded:\n")
print(donor_meta[, c("donor_short", "cohort", "pool")])

# ── Fix 1: donor_colours — built robustly from donor_meta ─────────────────────
# Always 8 entries: 6 donors + doublet + unassigned
all_donor_labels <- c(donor_meta$donor_short, "doublet", "unassigned")
colour_palette   <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                      "#59A14F", "#EDC948", "#FF0000", "#000000")

# Extend palette if more than 8 donors
if (length(all_donor_labels) > length(colour_palette)) {
  colour_palette <- c(colour_palette,
                      rainbow(length(all_donor_labels) - length(colour_palette)))
}

donor_colours <- setNames(
  colour_palette[seq_along(all_donor_labels)],
  all_donor_labels
)

cat("✔ Donor colours defined for:", paste(names(donor_colours), collapse = ", "), "\n")

# ── Fix 2: msigdbr column name — works with both old and new versions ──────────
get_ribo_genes <- function(c2_sets) {
  # Try new column name first, fall back to old
  gene_col <- if ("gene_symbol" %in% colnames(c2_sets)) {
    "gene_symbol"
  } else if ("human_gene_symbol" %in% colnames(c2_sets)) {
    "human_gene_symbol"
  } else {
    # Show available columns to help debug
    cat("Available msigdbr columns:", paste(colnames(c2_sets), collapse = ", "), "\n")
    stop("Cannot find gene symbol column in msigdbr output")
  }
  cat(glue("✔ Using msigdbr column: {gene_col}\n"))
  c2_sets[c2_sets$gs_name == "KEGG_RIBOSOME", gene_col, drop = TRUE]
}

# Compute c2_sets once
c2_sets  <- msigdbr(species = "Homo sapiens", category = "C2")
ribo_kegg <- get_ribo_genes(c2_sets)

# ══════════════════════════════════════════════════════════════════════════════
# LOOP OVER POOLS
# ══════════════════════════════════════════════════════════════════════════════
for (i in c("A", "B")) {
  
  cat(glue("\n{strrep('═', 60)}\n  QC Plots — Pool {i}\n{strrep('═', 60)}\n"))
  
  # ── Load raw CellRanger SCE ────────────────────────────────────────────────
  raw_file <- here("data", "SCEs", glue("Pool_{i}.CellRanger.SCE.rds"))
  sce_raw  <- if (file.exists(raw_file)) {
    sce_r <- readRDS(raw_file)
    sce_r$capture <- factor(sce_r$Sample)
    cat(glue("  Raw cells: {ncol(sce_r)}\n"))
    sce_r
  } else {
    cat(glue("  ⚠ Raw SCE not found — skipping barcode rank plot\n"))
    NULL
  }
  
  # ── Load preprocessed SCE ─────────────────────────────────────────────────
  pre_file <- here("data", "SCEs", glue("Pool_{i}.preprocessed.SCE.rds"))
  if (!file.exists(pre_file)) {
    cat(glue("  ⚠ Preprocessed SCE not found — run 03_preprocess.R first\n"))
    next
  }
  
  sce <- readRDS(pre_file)
  cat(glue("  Preprocessed cells: {ncol(sce)}\n"))
  
  # ── Recompute gene sets ────────────────────────────────────────────────────
  mito_set <- rownames(sce)[which(rowData(sce)$CHR == "MT")]
  ribo_set <- union(
    grep("^RP(S|L)", rownames(sce), value = TRUE),
    ribo_kegg
  )
  
  # Recompute QC metrics if missing
  if (!"subsets_Mito_percent" %in% colnames(colData(sce))) {
    cat("  Recomputing QC metrics...\n")
    sce <- addPerCellQC(
      sce,
      subsets = list(
        Mito = which(rownames(sce) %in% mito_set),
        Ribo = which(rownames(sce) %in% ribo_set)
      )
    )
  }
  
  # Donors present in this pool — only use colours for these
  # donors_in_pool <- intersect(
  #   levels(droplevels(factor(sce$genetic_donor))),
  #   names(donor_colours)
  # )
  donors_in_pool <- intersect(
    unique(as.character(sce$genetic_donor)),
    names(donor_colours)
  )
  pool_colours <- donor_colours[donors_in_pool]
  cat(glue("  Donors in pool: {paste(donors_in_pool, collapse=', ')}\n"))
  
  # ════════════════════════════════════════════════════════════════════════════
  # PLOT 1: Barcode rank plots
  # ════════════════════════════════════════════════════════════════════════════
  if (!is.null(sce_raw)) {
    pdf(here("output", "figures", "QC", glue("Pool_{i}_barcode_rank.pdf")),
        width = 12, height = 6)
    
    capture_names <- levels(sce_raw$capture)
    n_cap <- length(capture_names)
    par(mfrow = c(ceiling(n_cap / 2), 2), mar = c(4, 4, 3, 1))
    
    for (cn in capture_names) {
      sce_sub <- sce_raw[, sce_raw$capture == cn]
      tryCatch({
        bcrank <- barcodeRanks(counts(sce_sub))
        uniq   <- !duplicated(bcrank$rank)
        plot(
          x = bcrank$rank[uniq], y = bcrank$total[uniq],
          log = "xy", xlab = "Rank", ylab = "Total UMI",
          main = cn, pch = 16, cex = 0.3, col = "grey60"
        )
        abline(h = metadata(bcrank)$inflection,
               col = "darkgreen", lty = 2, lwd = 1.5)
        abline(h = metadata(bcrank)$knee,
               col = "dodgerblue", lty = 2, lwd = 1.5)
        legend("bottomleft",
               legend = c("Inflection", "Knee"),
               col    = c("darkgreen", "dodgerblue"),
               lty = 2, lwd = 1.5, cex = 0.8, bty = "n")
      }, error = function(e) {
        plot.new()
        title(main = paste(cn, "— error:", e$message))
      })
    }
    dev.off()
    cat("  ✔ barcode_rank.pdf\n")
  }
  
  # ════════════════════════════════════════════════════════════════════════════
  # PLOT 2: QC metrics violin per donor
  # ════════════════════════════════════════════════════════════════════════════
  meta_df <- as.data.frame(colData(sce)) %>%
    filter(!genetic_donor %in% c("doublet", "unassigned")) %>%
    mutate(
      genetic_donor = factor(
        as.character(genetic_donor),              # ← strip old factor levels first
        levels = unique(c(donor_meta$donor_short,  # ← unique() prevents duplicates
                          "doublet", "unassigned"))
      )
    )
  
  make_violin <- function(df, y_var, y_lab, title, colours,
                          hline = NULL) {
    p <- ggplot(df, aes_string(x = "genetic_donor", y = y_var,
                               fill = "genetic_donor")) +
      geom_violin(alpha = 0.7, scale = "width") +
      geom_boxplot(width = 0.12, outlier.size = 0.3,
                   fill = "white", alpha = 0.8) +
      scale_fill_manual(values = colours, guide = "none") +
      labs(y = y_lab, x = NULL, title = title) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title  = element_text(face = "bold"))
    if (!is.null(hline)) {
      p <- p + geom_hline(yintercept = hline,
                          linetype = "dashed", color = "firebrick",
                          linewidth = 0.8)
    }
    p
  }
  
  pdf(here("output", "figures", "QC", glue("Pool_{i}_QC_violin.pdf")),
      width = 14, height = 10)
  
  pv1 <- make_violin(meta_df %>% mutate(log_sum = log10(sum + 1)),
                     "log_sum", "log10(UMI + 1)", "Total UMI", pool_colours)
  pv2 <- make_violin(meta_df, "detected", "Genes detected",
                     "Genes detected", pool_colours)
  pv3 <- make_violin(meta_df, "subsets_Mito_percent", "Mito %",
                     "Mitochondrial %", pool_colours, hline = 20)
  pv4 <- make_violin(meta_df, "subsets_Ribo_percent", "Ribo %",
                     "Ribosomal %", pool_colours)
  
  print(pv1 + pv2 + pv3 + pv4 +
          plot_layout(ncol = 2) +
          plot_annotation(
            title    = glue("Pool {i}: QC metric distributions by donor"),
            subtitle = glue("{nrow(meta_df)} cells (doublets/unassigned excluded)"),
            theme    = theme(
              plot.title    = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, color = "grey40")
            )
          ))
  dev.off()
  cat("  ✔ QC_violin.pdf\n")
  
  # ════════════════════════════════════════════════════════════════════════════
  # PLOT 3: Mito% vs total counts scatter
  # ════════════════════════════════════════════════════════════════════════════
  pdf(here("output", "figures", "QC", glue("Pool_{i}_mito_scatter.pdf")),
      width = 9, height = 6)
  
  scatter_df <- as.data.frame(colData(sce)) %>%
    filter(!genetic_donor %in% c("doublet", "unassigned"))
  
  ps <- ggplot(scatter_df,
               aes(x = sum, y = subsets_Mito_percent,
                   color = genetic_donor)) +
    geom_point(size = 0.4, alpha = 0.5) +
    scale_x_log10() +
    scale_color_manual(values = pool_colours, name = "Donor") +
    geom_hline(yintercept = 20, linetype = "dashed",
               color = "firebrick", linewidth = 0.8) +
    labs(
      x     = "Total UMI (log10)",
      y     = "Mitochondrial %",
      title = glue("Pool {i}: Mito% vs total UMI"),
      subtitle = "Dashed line = 20% threshold"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
  
  print(ps)
  dev.off()
  cat("  ✔ mito_scatter.pdf\n")
  
  # ════════════════════════════════════════════════════════════════════════════
  # PLOT 4: Cells per donor bar chart
  # ════════════════════════════════════════════════════════════════════════════
  pdf(here("output", "figures", "QC", glue("Pool_{i}_donor_cells.pdf")),
      width = 9, height = 5)
  
  donor_counts <- as.data.frame(colData(sce)) %>%
    count(genetic_donor) %>%
    mutate(
      genetic_donor = factor(
        as.character(genetic_donor),              # ← strip old factor levels first
        levels = unique(c(donor_meta$donor_short,  # ← unique() prevents duplicates
                          "doublet", "unassigned"))
      )
    )  
  pd <- ggplot(donor_counts,
               aes(x = genetic_donor, y = n, fill = genetic_donor)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = donor_colours, guide = "none") +
    labs(
      x     = "Donor",
      y     = "Number of cells",
      title = glue("Pool {i}: Cells per donor after QC"),
      subtitle = glue("Total: {ncol(sce)} cells")
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          axis.text.x   = element_text(angle = 30, hjust = 1))
  
  print(pd)
  dev.off()
  cat("  ✔ donor_cells.pdf\n")
  
  # ════════════════════════════════════════════════════════════════════════════
  # PLOT 5: Doublet summary
  # ════════════════════════════════════════════════════════════════════════════
  if ("scDblFinder.class" %in% colnames(colData(sce))) {
    
    pdf(here("output", "figures", "QC", glue("Pool_{i}_doublets.pdf")),
        width = 9, height = 5)
    
    dbl_df <- as.data.frame(colData(sce)) %>%
      count(genetic_donor, scDblFinder.class) %>%
      group_by(genetic_donor) %>%
      mutate(pct = round(100 * n / sum(n), 1)) %>%
      ungroup()
    
    pdbl <- ggplot(dbl_df,
                   aes(x = genetic_donor, y = n,
                       fill = scDblFinder.class)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      geom_text(aes(label = paste0(pct, "%")),
                position = position_stack(vjust = 0.5),
                size = 3, color = "white", fontface = "bold") +
      scale_fill_manual(
        values = c("singlet" = "#2166AC", "doublet" = "#D73027"),
        name   = "Class"
      ) +
      labs(
        x     = "Donor",
        y     = "Cells",
        title = glue("Pool {i}: Singlet / doublet breakdown")
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title  = element_text(hjust = 0.5, face = "bold"),
            axis.text.x = element_text(angle = 30, hjust = 1))
    
    print(pdbl)
    dev.off()
    cat("  ✔ doublets.pdf\n")
  }
  
  # ── Summary stats ──────────────────────────────────────────────────────────
  cat(glue("\n  Pool {i} summary:\n"))
  cat(glue("  Total cells     : {ncol(sce)}\n"))
  cat(glue("  Median UMI      : {median(sce$sum)}\n"))
  cat(glue("  Median genes    : {median(sce$detected)}\n"))
  cat(glue("  Median mito%%   : {round(median(sce$subsets_Mito_percent), 1)}\n"))
  cat("  Cells per donor:\n")
  print(table(sce$genetic_donor))
}

cat(glue("\n{strrep('═', 60)}\n  ALL QC PLOTS COMPLETE\n{strrep('═', 60)}\n"))
cat(glue("  Output: {here('output', 'figures', 'QC')}\n\n"))
cat("  Files per pool:\n")
cat("    Pool_X_barcode_rank.pdf   — knee/inflection plots\n")
cat("    Pool_X_QC_violin.pdf      — violin distributions per donor\n")
cat("    Pool_X_mito_scatter.pdf   — mito% vs UMI scatter\n")
cat("    Pool_X_donor_cells.pdf    — cells per donor bar chart\n")
cat("    Pool_X_doublets.pdf       — singlet/doublet breakdown\n")

