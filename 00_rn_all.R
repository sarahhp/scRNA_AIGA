# ══════════════════════════════════════════════════════════════════════════════
# 00_run_all.R — Master pipeline runner
# ══════════════════════════════════════════════════════════════════════════════

setwd("/g/data/pq84/single_cell/10X/scRNA_auto/")

dyn.load("/home/577/mk9000/.local/lib/libglpk.so.40")
dyn.load("/g/data/pq84/mano/Rlib_v4.4.2/igraph/libs/igraph.so")
.libPaths("/g/data/pq84/mano/Rlib_v4.4.2")


scripts <- c(
  "scripts/01_load_cellranger.R",
  "scripts/02_demultiplex.R",
  "scripts/03_preprocess.R",
  "scripts/04_integrate_annotate.R",
  "scripts/05_DE_analysis.R"
)

for (script in scripts) {
  cat("\n══════════════════════════════════════════════\n")
  cat("Running:", script, "\n")
  cat("Time:", format(Sys.time()), "\n")
  cat("══════════════════════════════════════════════\n")
  source(script)
}

cat("\n✔ Pipeline complete:", format(Sys.time()), "\n")
