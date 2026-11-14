# What the pruning does to a panel at matched size, on the eleven TCGA cohorts. Same
# coding-only feature space, five outer folds, seed and ridge Cox learner as
# cardinality_curve.R, so the survfcbf, variance, CAR and MIM rows reproduce that curve
# exactly (paper/make_statistics.R checks it). Four questions:
#   1. redundancy within the panel: mean absolute Pearson correlation between the selected
#      genes, and the share of pairs with |r| >= 0.5, on the training rows of each fold;
#   2. mRMR at K = 5, 10 and 20 next to survfcbf top-K and the three rankers;
#   3. Cox on censoring: the panel is chosen on the true outcome, then a Cox model on the
#      panel is trained and scored with the status reversed (censoring as the event). Its
#      Harrell's C says how well the panel predicts follow-up length. A random ten-gene
#      panel, drawn from the fold's seeded stream, is the reference;
#   4. top-K at q = 4 and q = 10 (K = 10) against the default q = 2.
# It also runs MIM and mRMR at q = 2, the binning of top-K, so that the discretisation is
# not confounded with the pruning (by default they take q = 10 from Bommert's rule); optionally (RD_THETAS) stability at theta = 0.5 and 0.7
# besides 0.9 for the K = 10 panels; the selected panels per fold (<cohort>_panels.csv),
# from which those thresholds can be scored apart; and 100 random ten-gene
# panels per fold for the Cox on censoring (<cohort>_reverse_random.csv), so the real
# panel is placed in a distribution rather than against a single draw; and Uno's C
# truncated at the 80th percentile of each training fold's follow-up (uno_c_trunc, with
# the per-fold values in <cohort>_folds.csv, realdata/uno_truncated.R). The rows of the
# first run (selectors listed first, same folds and seeds) reproduce exactly.
#
# SLURM array, one task per cohort (cesga/submit_script.sbatch). The published run used
# the survfcbf build with the Spearman tie-break (CESGA, 2026-11-19, 35 to 60 min per
# cohort at 2 CPUs and 32 GB):
#   sbatch -p short --qos=short --time=03:00:00 --mem=32G --cpus-per-task=2 --array=1-11 \
#     cesga/submit_script.sbatch realdata/panel_redundancy.R
# Writes <cohort>.csv and <cohort>_reverse.csv to RD_OUTDIR (default
# results/realdata_redundancy). Inference: realdata/redundancy_inference.R.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners)
  library(stabm); library(farff); library(jsonlite); library(curl)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
stopifnot("this needs the tie-break build of survfcbf" = "tiebreak" %in% names(formals(survfcbf_select)))
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(getwd(), "realdata/Data"))
source("realdata/uno_truncated.R")
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
co <- names(BOMMERT_IDS)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_redundancy")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".csv"))
if (file.exists(out)) { cat(co, "exists, skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
time <- d$time; status <- d$status
cat(sprintf("%s: n=%d events=%d p=%d\n", co, length(time), sum(status), ncol(X)))
t0 <- proc.time()[["elapsed"]]

## ---- selectors, memoised per training fold (folds are shared across selectors) ----------
Ks <- c(5L, 10L, 20L)
cache <- new.env(parent = emptyenv())
memo <- function(key, ti, Xs, fn) {
  k <- sprintf("%s_%d_%.4f", key, nrow(Xs), sum(ti))
  if (is.null(cache[[k]])) assign(k, fn(), envir = cache)
  cache[[k]]
}
rank_of <- function(base, ti, st, Xs) memo(base, ti, Xs, function() switch(base,
  variance = apply(Xs, 2L, stats::var),
  carss    = abs(run_baseline(ent("carss"),    ti, st, Xs, prop = NULL, seed = 1L)$score),
  mim_surv = run_baseline(ent("mim_surv"), ti, st, Xs, prop = NULL, seed = 1L)$score))
topk_of <- function(q, ti, st, Xs) memo(paste0("fcbf_q", q), ti, Xs, function()
  survfcbf_select(ti, st, Xs, q = q, top = max(Ks))$selected)
mrmr_of <- function(ti, st, Xs) memo("mrmr", ti, Xs, function()     # greedy: the first K are a prefix
  run_baseline(ent("mrmr_surv"), ti, st, Xs, prop = max(Ks) / ncol(Xs), seed = 1L)$selected)

sels <- list()
for (K in Ks) local({ k <- K
  sels[[sprintf("survfcbf_q2_K%d", k)]] <<- function(ti, st, Xs) utils::head(topk_of(2L, ti, st, Xs), k)
  sels[[sprintf("mrmr_K%d", k)]]        <<- function(ti, st, Xs) utils::head(mrmr_of(ti, st, Xs), k)
  for (b in c("variance", "carss", "mim_surv")) local({ bb <- b
    sels[[sprintf("%s_K%d", bb, k)]] <<- function(ti, st, Xs) cut_ranking(rank_of(bb, ti, st, Xs), k = k) })
})
for (q in c(4L, 10L)) local({ qq <- q
  sels[[sprintf("survfcbf_q%d_K10", qq)]] <<- function(ti, st, Xs) utils::head(topk_of(qq, ti, st, Xs), 10L) })
## MIM and mRMR at the binning of top-K (q = 2), appended so the earlier selectors keep
## their order
mim_q2_of <- function(ti, st, Xs) memo("mim_q2", ti, Xs, function()
  mim_surv_select(ti, st, Xs, q = 2L)$score)
mrmr_q2_of <- function(ti, st, Xs) memo("mrmr_q2", ti, Xs, function()
  mrmr_surv_select(ti, st, Xs, prop = max(Ks) / ncol(Xs), q = 2L)$selected)
for (K in Ks) local({ k <- K
  sels[[sprintf("mim_surv_q2_K%d", k)]] <<- function(ti, st, Xs) cut_ranking(mim_q2_of(ti, st, Xs), k = k)
  sels[[sprintf("mrmr_q2_K%d", k)]]     <<- function(ti, st, Xs) utils::head(mrmr_q2_of(ti, st, Xs), k)
})
## extra stability thresholds for the K = 10 panels, off by default: at theta = 0.5 the
## chance correction of stabm on ~16k coding genes takes hours per panel (the first run with
## them timed out at 4 h in every cohort), so they run apart, from the saved panels
THETAS <- as.numeric(strsplit(Sys.getenv("RD_THETAS", ""), ",", fixed = TRUE)[[1]])

## in-panel redundancy, on the training rows of each fold
redundancy <- function(bm) {
  v <- vapply(bm$folds, function(fd) {
    f <- as.character(fd$features)
    if (length(f) < 2L) return(c(NA_real_, NA_real_))
    r <- abs(stats::cor(X[fd$train_set, f, drop = FALSE]))
    u <- r[upper.tri(r)]
    c(mean(u), mean(u >= 0.5))
  }, numeric(2))
  c(mean_abs_r = mean(v[1, ], na.rm = TRUE), share_r05 = mean(v[2, ], na.rm = TRUE))
}

L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")
task <- surv_task(time, status, X)
sim_mat <- NULL
panels <- list(); fold_rows <- list()
rows <- lapply(names(sels), function(nm) {
  bm <- nested_cv_select(task, sels[[nm]], L$clone(deep = TRUE), outer_folds = 5L, seed = 1L)
  sc <- score_benchmark(bm)
  if (is.null(sim_mat)) sim_mat <<- similarity_matrix(bm, scope = "full")
  st <- benchmark_stability(bm, theta = 0.9, seed = 1L, sim_mat = sim_mat)
  st_extra <- vapply(c(0.5, 0.7), function(th) if (th %in% THETAS && grepl("_K10$", nm))
    benchmark_stability(bm, theta = th, seed = 1L, sim_mat = sim_mat)$stability else NA_real_, 0)
  red <- redundancy(bm)
  uf <- uno_by_fold(bm)
  fold_rows[[nm]] <<- data.frame(cohort = co, selector = nm, uf)
  panels[[nm]] <<- do.call(rbind, lapply(seq_along(bm$folds), function(f) {
    ft <- as.character(bm$folds[[f]]$features)
    if (!length(ft)) return(NULL)
    data.frame(cohort = co, selector = nm, fold = f, position = seq_along(ft), feature = ft)
  }))
  s <- stats::setNames(sc$summary$mean, sc$summary$metric)
  data.frame(cohort = co, n = length(time), events = sum(status), p = ncol(X), selector = nm,
             uno_c = s[["uno_c"]], harrell_c = s[["harrell_c"]], ibs = s[["ibs"]],
             uno_c_se = sc$summary$se[sc$summary$metric == "uno_c"],
             stability = st$stability, parsimony = st$parsimony_mean,
             mean_abs_r = red[["mean_abs_r"]], share_r05 = red[["share_r05"]],
             n_failed = sc$n_failed, n_null_fit = sc$n_null_fit,
             stability_t05 = st_extra[1], stability_t07 = st_extra[2],
             uno_c_trunc = mean(uf$uno_c_trunc), row.names = NULL)
})
res <- do.call(rbind, rows)

## ---- Cox on censoring: panel chosen on the true outcome, learner trained and scored on the
## reversed one
task_rev <- surv_task(time, 1L - status, X)
REV <- c("survfcbf_q2_K10", "survfcbf_q4_K10", "survfcbf_q10_K10", "variance_K10", "carss_K10",
         "mim_surv_K10", "mrmr_K10", "mim_surv_q2_K10", "mrmr_q2_K10")
rev_sels <- lapply(sels[REV], function(f) { force(f); function(ti, st, Xs) f(ti, 1L - st, Xs) })
rev_sels$random_K10 <- function(ti, st, Xs) sample.int(ncol(Xs), 10L)   # seeded by the fold stream
rev <- do.call(rbind, lapply(names(rev_sels), function(nm) {
  bm <- nested_cv_select(task_rev, rev_sels[[nm]], L$clone(deep = TRUE), outer_folds = 5L, seed = 1L)
  sc <- score_benchmark(bm)
  data.frame(cohort = co, selector = nm,
             harrell_c_censoring = sc$summary$mean[sc$summary$metric == "harrell_c"],
             harrell_c_censoring_se = sc$summary$se[sc$summary$metric == "harrell_c"],
             n_failed = sc$n_failed, n_null_fit = sc$n_null_fit, row.names = NULL)
}))

## 100 random ten-gene panels per fold. Each draw has its own seed (replicate and training
## set), taken under with_seed() so the global RNG that the learner then uses is left as it
## was: the random panel is the only thing that changes between replicates.
R_RAND <- 100L
rand <- do.call(rbind, lapply(seq_len(R_RAND), function(r) {
  f <- function(ti, st, Xs) survfcbf:::with_seed(
    as.integer((r * 104729 + round(sum(ti))) %% 2147483647), sample.int(ncol(Xs), 10L))
  bm <- nested_cv_select(task_rev, f, L$clone(deep = TRUE), outer_folds = 5L, seed = 1L)
  sc <- score_benchmark(bm)
  data.frame(cohort = co, replicate = r,
             harrell_c_censoring = sc$summary$mean[sc$summary$metric == "harrell_c"],
             n_failed = sc$n_failed, n_null_fit = sc$n_null_fit, row.names = NULL)
}))

utils::write.csv(res, out, row.names = FALSE)
utils::write.csv(rev, file.path(outdir, paste0(co, "_reverse.csv")), row.names = FALSE)
utils::write.csv(rand, file.path(outdir, paste0(co, "_reverse_random.csv")), row.names = FALSE)
utils::write.csv(do.call(rbind, panels), file.path(outdir, paste0(co, "_panels.csv")), row.names = FALSE)
utils::write.csv(do.call(rbind, fold_rows), file.path(outdir, paste0(co, "_folds.csv")), row.names = FALSE)
## versions actually loaded, and when the survfcbf build was installed (NA under load_all)
key <- c("survfcbf", "glmnet", "mlr3", "mlr3proba", "mlr3extralearners", "carSurv", "praznik",
         "stabm", "survival")
built <- utils::packageDescription("survfcbf")$Built
utils::write.csv(data.frame(package = key,
                            version = vapply(key, function(k) as.character(utils::packageVersion(k)), ""),
                            R = R.version.string, survfcbf_built = if (is.null(built)) NA_character_ else built),
                 file.path(outdir, paste0("versions_", co, ".csv")), row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, paste0("sessionInfo_", co, ".txt")))
print(res[, c("selector", "uno_c", "stability", "parsimony", "mean_abs_r")], row.names = FALSE)
print(rev, row.names = FALSE)
cat(sprintf("%s done in %.1f min\nPANEL_REDUNDANCY_DONE\n", co, (proc.time()[["elapsed"]] - t0) / 60))
