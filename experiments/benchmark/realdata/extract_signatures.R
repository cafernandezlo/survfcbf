# Extract the FINAL (full-cohort) signatures for the enrichment analysis.
# SLURM array: task i handles cohort i, refits every
# selector on the WHOLE cohort (no CV -- the signature a practitioner would deploy) and
# records the selected feature names at matched cardinalities K in {5,10,20} (the
# sensitivity range). survfcbf_topk is the redesign (q=2, top-K); the rankers (variance,
# carss, mim) are cut to the same K for a fair, cardinality-matched overlap. One CSV +
# background per cohort, skip-if-exists, so a timed-out or partial run just re-runs the
# missing cohorts (the per-cohort carss fit at p~23k is ~10 min, hence one task each).
# Measured feature names are dumped so the ID format is verified data-driven. Run
# in the apptainer container on CESGA (survfcbf must carry the `top` arg).
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl) })
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R")
Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))

Ks      <- c(5L, 10L, 20L)
rankers <- c("variance", "carss", "mim_surv")
reg     <- baseline_registry()
ent     <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot("SLURM_ARRAY_TASK_ID out of range" = !is.na(id) && id >= 1L && id <= length(cohorts))
co <- cohorts[id]

coding_only <- nzchar(Sys.getenv("CODING_ONLY", ""))     # protein-coding-restricted sub-analysis
outdir <- file.path("results/realdata_signatures", if (coding_only) "cells_coding" else "cells")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(outdir, paste0(co, ".csv"))
bgf <- file.path(outdir, paste0(co, "_bg.rds"))
if (file.exists(csv) && file.exists(bgf)) { cat(co, "exists, skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
if (coding_only) {                                       # keep only protein-coding rna features (HGNC)
  pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
  keep <- grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc
  X <- X[, keep, drop = FALSE]; nm <- colnames(X)
  cat(sprintf("  coding-only: %d -> %d features\n", ncol(d$X), ncol(X)))
}
saveRDS(nm, bgf)
sig <- list()
add <- function(sel, K, idx)
  sig[[length(sig) + 1L]] <<- data.frame(cohort = co, selector = sel, K = K,
    rank = seq_along(idx), feature = nm[idx], stringsAsFactors = FALSE)
topk20 <- survfcbf_select(d$time, d$status, X, q = 2L, top = max(Ks))$selected  # top-K is a prefix
for (K in Ks) add("survfcbf_topk", K, head(topk20, K))
for (rk in rankers) {
  sc <- run_baseline(ent(rk), d$time, d$status, X, prop = NULL, seed = 1L)$score
  for (K in Ks) add(rk, K, cut_ranking(sc, k = K))
}
utils::write.csv(do.call(rbind, sig), csv, row.names = FALSE)
cat(sprintf("%s done: n=%d ev=%d p=%d\n  topk@10: %s\n", co, length(d$time), sum(d$status),
            ncol(X), paste(nm[head(topk20, 10L)], collapse = " ")))
cat("CELL_DONE\n")
