# Repeated cross-validation of the METABRIC comparison: top-K (q = 2, ten genes) against
# the three rankers at the ten-percent cut-off, on the complete expression profile, as in
# metabric_benchmark.R. Ten five-fold splits (seeds 1 to 10, seed 1 being the published
# split, which it must reproduce), scoring Uno's C untruncated and truncated at the 80th
# percentile of each training fold's follow-up (realdata/uno_truncated.R), the integrated
# Brier score and stability. The design matches repeated_cv_main.R on TCGA.
#
# Runs locally (METABRIC_DIR as in load_metabric.R). RCV_REPS lists the repetitions this
# process runs (default 1:10), so several processes can share them; each writes rep<r>.csv
# to results/realdata_metabric/repeated_cv.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners); library(stabm)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_metabric.R"); source("realdata/uno_truncated.R")
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
OUT <- "results/realdata_metabric/repeated_cv"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
reps <- eval(parse(text = Sys.getenv("RCV_REPS", "1:10")))
PROP <- 0.1; TOPK <- 10L                                 # as cfg in metabric_benchmark.R
quiet <- function(expr) withCallingHandlers(expr, survfcbf_followup_warning = function(w) invokeRestart("muffleWarning"))
sels <- list(
  survfcbf_topk = function(ti, st, X) quiet(survfcbf_select(ti, st, X, q = 2L, top = TOPK)$selected),
  variance = function(ti, st, X) run_baseline(ent("variance"), ti, st, X, prop = PROP)$selected,
  carss    = function(ti, st, X) run_baseline(ent("carss"), ti, st, X, prop = PROP)$selected,
  mim_surv = function(ti, st, X) run_baseline(ent("mim_surv"), ti, st, X, prop = PROP)$selected)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")
d <- load_metabric(); task <- surv_task(d$time, d$status, d$X)
sim_mat <- NULL
for (r in reps) {
  out <- file.path(OUT, sprintf("rep%02d.csv", r))
  if (file.exists(out)) { cat(basename(out), "exists, skip\n"); next }
  t0 <- proc.time()[["elapsed"]]
  res <- do.call(rbind, lapply(names(sels), function(nm) {
    bm <- nested_cv_select(task, sels[[nm]], L$clone(deep = TRUE), outer_folds = 5L, seed = r)
    sc <- score_benchmark(bm)
    if (is.null(sim_mat)) sim_mat <<- similarity_matrix(bm, scope = "full")
    st <- benchmark_stability(bm, theta = 0.9, seed = r, sim_mat = sim_mat)
    s <- stats::setNames(sc$summary$mean, sc$summary$metric)
    data.frame(cohort = "METABRIC", repetition = r, seed = r, selector = nm,
               uno_c = s[["uno_c"]], uno_c_trunc = mean(uno_by_fold(bm)$uno_c_trunc),
               harrell_c = s[["harrell_c"]], ibs = s[["ibs"]], stability = st$stability,
               parsimony = st$parsimony_mean, n_failed = sc$n_failed, n_null_fit = sc$n_null_fit,
               row.names = NULL)
  }))
  utils::write.csv(res, out, row.names = FALSE)
  cat(sprintf("METABRIC rep %d done in %.1f min\n", r, (proc.time()[["elapsed"]] - t0) / 60))
}
cat("METABRIC_REPEATED_CV_DONE\n")
