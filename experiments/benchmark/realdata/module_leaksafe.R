# Leakage-safe module-FCBF stability (module pivot finalisation). The
# quick-check probes fixed the modules on the full cohort, which makes selection look
# artificially stable. Here modules are REBUILT on each subsample (fully leakage-safe), so
# this is the honest stability number. The selected units are the genes in
# the top-K modules (union), compared across subsamples with the Nogueira index -- i.e.
# "do the same genes get selected (via modules) across resamples". gene-level survfcbf and
# variance are already leakage-safe and re-run on the same subsamples for a direct contrast.
# NOTE: the definitive module-selection metric (cross-subsample module matching by gene
# overlap) is a method-writing decision; the gene-union Nogueira here is the pragmatic,
# defensible leakage-safe read. SLURM array, coding-only.
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl); library(stabm) })
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")); co <- cohorts[id]
outdir <- "results/realdata_modules_leaksafe"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".rds")); if (file.exists(out)) { cat(co, "skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
cat(sprintf("%s: n=%d ev=%d p=%d\n", co, length(d$time), sum(d$status), ncol(X)))
K <- 10L; B <- 15L
su_target <- function(t, s) discretize_quantile(survival_target(t, s, method = "martingale"), q = 2L)$binned

geneset <- function(mi, rel) sort(unique(unlist(mi[order(rel, decreasing = TRUE)[seq_len(K)]])))

set.seed(1)
Lsu <- Lcar <- Lgene <- Lvar <- vector("list", B)
for (b in seq_len(B)) {
  i <- sample(length(d$time), round(0.8 * length(d$time)))
  Xs <- X[i, , drop = FALSE]; ts <- d$time[i]; ss <- d$status[i]
  hc  <- stats::hclust(stats::as.dist(1 - abs(stats::cor(Xs))), method = "complete")   # rebuilt per subsample
  mod <- stats::cutree(hc, k = 400L); big <- as.integer(names(which(table(mod) >= 5L)))
  mi  <- lapply(big, function(m) which(mod == m))
  EG  <- vapply(mi, function(g) prcomp(Xs[, g, drop = FALSE], scale. = FALSE)$x[, 1], numeric(nrow(Xs)))
  yb  <- su_target(ts, ss)
  rel_su  <- apply(EG, 2L, function(e) symmetrical_uncertainty(discretize_quantile(e, q = 2L)$binned, yb))
  rel_car <- abs(run_baseline(ent("carss"), ts, ss, EG, prop = NULL, seed = 1L)$score)
  Lsu[[b]]   <- geneset(mi, rel_su)
  Lcar[[b]]  <- geneset(mi, rel_car)
  Lgene[[b]] <- survfcbf_select(ts, ss, Xs, q = 2L, top = K)$selected
  Lvar[[b]]  <- cut_ranking(apply(Xs, 2L, stats::var), k = K)
  cat(sprintf("  subsample %d/%d\n", b, B))
}
res <- data.frame(cohort = co, n = length(d$time), events = sum(d$status), p = ncol(X),
  genes_per_sig = round(mean(lengths(Lsu)), 0),
  stab_module_su_ls  = stabm::stabilityNogueira(Lsu,   p = ncol(X)),
  stab_module_car_ls = stabm::stabilityNogueira(Lcar,  p = ncol(X)),
  stab_gene          = stabm::stabilityNogueira(Lgene, p = ncol(X)),
  stab_var           = stabm::stabilityNogueira(Lvar,  p = ncol(X)), stringsAsFactors = FALSE)
saveRDS(res, out)
cat(sprintf("%s: module-SU-ls %.3f | module-CAR-ls %.3f | gene %.3f | variance %.3f (sig ~%d genes)\n",
            co, res$stab_module_su_ls, res$stab_module_car_ls, res$stab_gene, res$stab_var, res$genes_per_sig))
cat("MODULE_LEAKSAFE_DONE\n")
