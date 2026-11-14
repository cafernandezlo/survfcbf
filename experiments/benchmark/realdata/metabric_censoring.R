# Cox model on censoring in METABRIC, the negative control for the follow-up check of
# panel_redundancy.R. METABRIC has long follow-up and events are more than half of the
# patients, so the upper class of the q = 2 martingale target holds few censored patients
# and a panel chosen on it should carry little follow-up information. Same design as in
# TCGA: each ten-gene panel is chosen on the true outcome, then the ridge Cox model is
# trained and scored with the status reversed (five folds, seed 1), against 100 random
# ten-gene panels per fold. Panels: top-K at q = 2 and q = 10, MIM and mRMR (praznik,
# their default q = 10), on the complete expression profile of metabric_benchmark.R.
#
# Runs locally: METABRIC_DIR must point at the folder of the two snapshots (see
# load_metabric.R). Writes censoring_cox.csv and censoring_cox_random.csv to
# results/realdata_metabric.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_metabric.R")
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
OUT <- "results/realdata_metabric"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

d <- load_metabric(); time <- d$time; status <- d$status; X <- d$X
ym <- survival_target(time, status, method = "martingale")
share_censored_upper <- mean(status[ym > stats::median(ym)] == 0L)
cat(sprintf("METABRIC: n=%d events=%d p=%d, censored share of the q = 2 upper class %.3f\n",
            length(time), sum(status), ncol(X), share_censored_upper))

K <- 10L
quiet <- function(expr) withCallingHandlers(expr, survfcbf_followup_warning = function(w) invokeRestart("muffleWarning"))
cache <- new.env(parent = emptyenv())
memo <- function(key, ti, Xs, fn) {
  k <- sprintf("%s_%d_%.4f", key, nrow(Xs), sum(ti))
  if (is.null(cache[[k]])) assign(k, fn(), envir = cache)
  cache[[k]]
}
sels <- list(
  survfcbf_q2_K10 = function(ti, st, Xs) memo("q2", ti, Xs, function() quiet(survfcbf_select(ti, st, Xs, q = 2L, top = K)$selected)),
  survfcbf_q10_K10 = function(ti, st, Xs) memo("q10", ti, Xs, function() quiet(survfcbf_select(ti, st, Xs, q = 10L, top = K)$selected)),
  mim_surv_K10 = function(ti, st, Xs) memo("mim", ti, Xs, function()
    cut_ranking(run_baseline(ent("mim_surv"), ti, st, Xs, prop = NULL, seed = 1L)$score, k = K)),
  mrmr_K10 = function(ti, st, Xs) memo("mrmr", ti, Xs, function()
    run_baseline(ent("mrmr_surv"), ti, st, Xs, prop = K / ncol(Xs), seed = 1L)$selected))
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")
task_rev <- surv_task(time, 1L - status, X)
score_rev <- function(f) {
  bm <- nested_cv_select(task_rev, f, L$clone(deep = TRUE), outer_folds = 5L, seed = 1L)
  sc <- score_benchmark(bm)
  c(harrell_c_censoring = sc$summary$mean[sc$summary$metric == "harrell_c"],
    harrell_c_censoring_se = sc$summary$se[sc$summary$metric == "harrell_c"],
    n_failed = sc$n_failed, n_null_fit = sc$n_null_fit)
}
res <- do.call(rbind, lapply(names(sels), function(nm) {
  f <- sels[[nm]]
  data.frame(selector = nm, t(score_rev(function(ti, st, Xs) f(ti, 1L - st, Xs))))
}))
rand <- do.call(rbind, lapply(seq_len(100L), function(r) {
  f <- function(ti, st, Xs) survfcbf:::with_seed(
    as.integer((r * 104729 + round(sum(ti))) %% 2147483647), sample.int(ncol(Xs), K))
  data.frame(replicate = r, t(score_rev(f)))
}))
res$random_mean <- mean(rand$harrell_c_censoring)
res$percentile_among_random <- vapply(res$harrell_c_censoring, function(v) mean(rand$harrell_c_censoring < v), 0)
res$share_censored_upper <- share_censored_upper
res$n <- length(time); res$events <- sum(status)
utils::write.csv(res, file.path(OUT, "censoring_cox.csv"), row.names = FALSE)
utils::write.csv(rand, file.path(OUT, "censoring_cox_random.csv"), row.names = FALSE)
print(res, row.names = FALSE)
cat(sprintf("random panels: mean %.3f, 2.5-97.5%% %.3f to %.3f\n", mean(rand$harrell_c_censoring),
            stats::quantile(rand$harrell_c_censoring, 0.025), stats::quantile(rand$harrell_c_censoring, 0.975)))
cat(sprintf("done in %.1f min\nMETABRIC_CENSORING_DONE\n", (proc.time()[["elapsed"]] - t0) / 60))
