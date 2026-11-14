# TEST 2: does PARSIMONY cost prediction? (the evidence gap). Sweeps the
# selected-set size K and traces downstream Uno's C, so the question "is the predictive
# signal in a sparse subset, or distributed across many genes?" gets a direct answer. A
# flat or rising curve means sparse selection discards information and the honest guidance
# for genomic survival is to use dense/regularised or aggregated models -- which would also
# explain why the ultra-parsimonious survfcbf cannot win on prediction by construction.
#
# Ranking selectors sweep the full K range; survfcbf_topk only to K=100 (its Markov-blanket
# pruning is O(p*k) and top-K is a prefix, so one call at max K is cut down). Per-fold
# rankings are memoised so carss is fitted once per fold, not once per K. Same harness,
# folds and measures as the main benchmark. SLURM array, coding-only.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners)
  library(stabm); library(farff); library(jsonlite); library(curl)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")); co <- cohorts[id]
## Partial reruns: RD_BASES (comma-separated: variance, carss, mim_surv, survfcbf) runs a
## subset and RD_OUTDIR writes elsewhere, so a rerun never overwrites the full run. A
## base's rows do not depend on the others: folds, seeds and similarity matrix are shared.
BASES  <- strsplit(Sys.getenv("RD_BASES", "variance,carss,mim_surv,survfcbf"), ",", fixed = TRUE)[[1]]
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_cardinality")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".csv")); if (file.exists(out)) { cat(co, "skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
cat(sprintf("%s: n=%d ev=%d p=%d\n", co, length(d$time), sum(d$status), ncol(X)))

Ks     <- c(5L, 10L, 20L, 50L, 100L, 500L, 2000L)
Ks_fcbf <- Ks[Ks <= 100L]
cache  <- new.env(parent = emptyenv())
memo <- function(key, ti, X, fn) {                     # folds are shared, so cache per fold
  k <- sprintf("%s_%d_%.4f", key, nrow(X), sum(ti))
  if (!is.null(cache[[k]])) return(cache[[k]])
  v <- fn(); assign(k, v, envir = cache); v
}
rank_of <- function(base, ti, st, Xs) memo(base, ti, Xs, function() switch(base,
  variance = apply(Xs, 2L, stats::var),
  carss    = abs(run_baseline(ent("carss"),    ti, st, Xs, prop = NULL, seed = 1L)$score),
  mim_surv = run_baseline(ent("mim_surv"), ti, st, Xs, prop = NULL, seed = 1L)$score))

sels <- list()
for (base in intersect(c("variance", "carss", "mim_surv"), BASES)) for (K in Ks)
  sels[[sprintf("%s_K%d", base, K)]] <- local({ b <- base; k <- K
    function(ti, st, Xs) cut_ranking(rank_of(b, ti, st, Xs), k = k) })
if ("survfcbf" %in% BASES) for (K in Ks_fcbf)          # top-K is a prefix: one call, then cut
  sels[[sprintf("survfcbf_K%d", K)]] <- local({ k <- K
    function(ti, st, Xs) head(memo("fcbf", ti, Xs, function()
      survfcbf_select(ti, st, Xs, q = 2L, top = max(Ks_fcbf))$selected), k) })
stopifnot("RD_BASES names no known base" = length(sels) > 0L)

task <- surv_task(d$time, d$status, X)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")   # 1se -> null model (?score_benchmark)
tab <- benchmark_selectors(task, sels, L$clone(deep = TRUE), outer_folds = 5L,
                           manifest = NULL, seed = 1L)
sp <- strsplit(tab$selector, "_K", fixed = TRUE)
res <- data.frame(cohort = co, n = length(d$time), events = sum(d$status), p = ncol(X),
                  base = vapply(sp, `[`, character(1), 1L),
                  K = as.integer(vapply(sp, `[`, character(1), 2L)),
                  tab[, setdiff(names(tab), "selector")], row.names = NULL, check.names = FALSE)
utils::write.csv(res, out, row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, sprintf("sessionInfo_%s.txt", co)))
w <- reshape(res[, c("base", "K", "uno_c")], idvar = "K", timevar = "base", direction = "wide")
cat("=== Uno's C vs K ===\n"); print(w[order(w$K), ], row.names = FALSE)
cat("CARDINALITY_DONE\n")
