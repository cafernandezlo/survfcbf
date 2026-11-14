# Sensitivity of top-K to the survival target. The default target, martingale residuals
# split at the median (q = 2), puts the censored patients with the shortest follow-up in
# the upper class when events are fewer than half the patients, so relevance can partly
# reflect follow-up. The IPCW alternative of survfcbf_select() (target_method =
# "ipcw_binned") bins the observed time and weights the patients by inverse probability of
# censoring, censored patients carrying zero weight, so follow-up does not enter that way.
# For each cohort, top-K (q = 2, K = 10) is run with both targets on the coding-only
# features and the five folds of panel_redundancy.R (seed 1): Uno's C untruncated and
# truncated, stability, redundancy within the panel, the overlap of the two panels per fold
# (Jaccard), and the Cox model on censoring of panel_redundancy.R (Harrell's C with the
# status reversed, the panel chosen on the true outcome).
#
# SLURM array, one task per cohort:
#   sbatch -p short --qos=short --time=03:00:00 --mem=32G --cpus-per-task=1 --array=1-11 \
#     cesga/submit_script.sbatch realdata/target_sensitivity.R
# Writes <cohort>.csv to RD_OUTDIR (default results/realdata_target_sensitivity).
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners)
  library(stabm); library(farff); library(jsonlite); library(curl)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
stopifnot("this needs the tie-break build of survfcbf" = "tiebreak" %in% names(formals(survfcbf_select)))
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(getwd(), "realdata/Data"))
source("realdata/uno_truncated.R")
co <- names(BOMMERT_IDS)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_target_sensitivity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".csv"))
if (file.exists(out)) { cat(co, "exists, skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
time <- d$time; status <- d$status
t0 <- proc.time()[["elapsed"]]
## upper class of the q = 2 martingale target on the whole cohort, for reference
ym <- survival_target(time, status, method = "martingale")
upper <- ym > stats::median(ym)
share_censored_upper <- mean(status[upper] == 0L)

K <- 10L
quiet <- function(expr) withCallingHandlers(expr, survfcbf_followup_warning = function(w) invokeRestart("muffleWarning"))
sels <- list(
  martingale = function(ti, st, Xs) quiet(survfcbf_select(ti, st, Xs, q = 2L, top = K)$selected),
  ipcw       = function(ti, st, Xs) survfcbf_select(ti, st, Xs, q = 2L, top = K, target_method = "ipcw_binned")$selected)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")
task <- surv_task(time, status, X); task_rev <- surv_task(time, 1L - status, X)
sim_mat <- NULL; panels <- list()
rows <- lapply(names(sels), function(nm) {
  bm <- nested_cv_select(task, sels[[nm]], L$clone(deep = TRUE), outer_folds = 5L, seed = 1L)
  sc <- score_benchmark(bm)
  if (is.null(sim_mat)) sim_mat <<- similarity_matrix(bm, scope = "full")
  st <- benchmark_stability(bm, theta = 0.9, seed = 1L, sim_mat = sim_mat)
  panels[[nm]] <<- lapply(bm$folds, function(fd) as.character(fd$features))
  red <- mean(vapply(bm$folds, function(fd) {
    f <- as.character(fd$features); if (length(f) < 2L) return(NA_real_)
    r <- abs(stats::cor(X[fd$train_set, f, drop = FALSE])); mean(r[upper.tri(r)]) }, 0), na.rm = TRUE)
  ## the same selector on the reversed task: the panel is chosen on the true status
  f <- sels[[nm]]; rsel <- function(ti, st, Xs) f(ti, 1L - st, Xs)
  bmr <- nested_cv_select(task_rev, rsel, L$clone(deep = TRUE), outer_folds = 5L, seed = 1L)
  scr <- score_benchmark(bmr)
  s <- stats::setNames(sc$summary$mean, sc$summary$metric)
  data.frame(cohort = co, n = length(time), events = sum(status), p = ncol(X), target = nm,
             share_censored_upper = share_censored_upper,
             uno_c = s[["uno_c"]], uno_c_trunc = mean(uno_by_fold(bm)$uno_c_trunc),
             harrell_c = s[["harrell_c"]], ibs = s[["ibs"]], stability = st$stability,
             parsimony = st$parsimony_mean, mean_abs_r = red,
             harrell_c_censoring = scr$summary$mean[scr$summary$metric == "harrell_c"],
             n_failed = sc$n_failed, n_null_fit = sc$n_null_fit,
             n_null_fit_censoring = scr$n_null_fit, row.names = NULL)
})
res <- do.call(rbind, rows)
res$jaccard_between_targets <- mean(mapply(function(a, b) length(intersect(a, b)) / length(union(a, b)),
                                           panels$martingale, panels$ipcw))
utils::write.csv(res, out, row.names = FALSE)
print(res, row.names = FALSE)
cat(sprintf("%s done in %.1f min\nTARGET_SENSITIVITY_DONE\n", co, (proc.time()[["elapsed"]] - t0) / 60))
