# Repeated cross-validation of the matched-size comparisons. The published runs use one
# five-fold split per cohort (seed 1), so a paired difference between two selectors carries
# the resampling noise of that single split. This repeats the split ten times (seeds 1 to
# 10, seed 1 being the published split) for top-K at q = 2 and its rivals at K genes (RCV_K,
# default 10): mRMR and MIM at their default q = 10 and at q = 2, CAR and the variance
# filter. Same coding-only features, harness, learner and measures as panel_redundancy.R.
#
# SLURM array over (cohort, repetition) units, with the repetitions of a cohort kept
# together so the data and the similarity matrix load once per task:
#   sbatch -p short --qos=short --time=03:00:00 --mem=32G --cpus-per-task=2 --array=1-55 \
#     --export=ALL,N_TASKS=55 cesga/submit_script.sbatch realdata/repeated_cv.R
# Add RCV_K=5 or RCV_K=20 to --export for the other panel sizes. Writes <cohort>_rep<r>.csv
# to RD_OUTDIR (default results/realdata_repeated_cv for K = 10, otherwise
# results/realdata_repeated_cv_K<K>).
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
K <- as.integer(Sys.getenv("RCV_K", "10"))
stopifnot("RCV_K must be a positive integer" = !is.na(K) && K > 0L)
outdir <- Sys.getenv("RD_OUTDIR", paste0("results/realdata_repeated_cv", if (K == 10L) "" else paste0("_K", K)))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

REPS <- as.integer(Sys.getenv("RCV_REPS", "10"))
units <- expand.grid(rep = seq_len(REPS), cohort = names(BOMMERT_IDS), stringsAsFactors = FALSE)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_tasks <- as.integer(Sys.getenv("N_TASKS", as.character(nrow(units))))
per <- ceiling(nrow(units) / n_tasks)
mine <- seq.int((id - 1L) * per + 1L, min(id * per, nrow(units)))
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")

loaded <- ""
for (i in mine) {
  co <- units$cohort[i]; r <- units$rep[i]
  out <- file.path(outdir, sprintf("%s_rep%02d.csv", co, r))
  if (file.exists(out)) { cat(basename(out), "exists, skip\n"); next }
  t0 <- proc.time()[["elapsed"]]
  if (loaded != co) {
    d <- load_bommert(co); X <- d$X; nm <- colnames(X)
    pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
    X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
    task <- surv_task(d$time, d$status, X)
    sim_mat <- NULL                         # outcome-free and fold-free: one per cohort
    loaded <- co
  }
  cache <- new.env(parent = emptyenv())     # folds differ between repetitions
  memo <- function(key, ti, Xs, fn) {
    k <- sprintf("%s_%d_%.4f", key, nrow(Xs), sum(ti))
    if (is.null(cache[[k]])) assign(k, fn(), envir = cache)
    cache[[k]]
  }
  sels <- list(
    survfcbf_q2 = function(ti, st, Xs) utils::head(memo("fcbf", ti, Xs, function()
      survfcbf_select(ti, st, Xs, q = 2L, top = K)$selected), K),
    mrmr = function(ti, st, Xs) run_baseline(ent("mrmr_surv"), ti, st, Xs, prop = K / ncol(Xs), seed = 1L)$selected,
    mrmr_q2 = function(ti, st, Xs) mrmr_surv_select(ti, st, Xs, prop = K / ncol(Xs), q = 2L)$selected,
    mim_surv = function(ti, st, Xs) cut_ranking(run_baseline(ent("mim_surv"), ti, st, Xs, prop = NULL, seed = 1L)$score, k = K),
    mim_surv_q2 = function(ti, st, Xs) cut_ranking(mim_surv_select(ti, st, Xs, q = 2L)$score, k = K),
    carss = function(ti, st, Xs) cut_ranking(abs(run_baseline(ent("carss"), ti, st, Xs, prop = NULL, seed = 1L)$score), k = K),
    variance = function(ti, st, Xs) cut_ranking(apply(Xs, 2L, stats::var), k = K))
  res <- do.call(rbind, lapply(names(sels), function(nm) {
    bm <- nested_cv_select(task, sels[[nm]], L$clone(deep = TRUE), outer_folds = 5L, seed = r)
    sc <- score_benchmark(bm)
    if (is.null(sim_mat)) sim_mat <<- similarity_matrix(bm, scope = "full")
    st <- benchmark_stability(bm, theta = 0.9, seed = r, sim_mat = sim_mat)
    s <- stats::setNames(sc$summary$mean, sc$summary$metric)
    uf <- uno_by_fold(bm)
    data.frame(cohort = co, repetition = r, seed = r, selector = paste0(nm, "_K", K),
               uno_c = s[["uno_c"]], uno_c_trunc = mean(uf$uno_c_trunc),
               harrell_c = s[["harrell_c"]], ibs = s[["ibs"]],
               stability = st$stability, parsimony = st$parsimony_mean,
               n_failed = sc$n_failed, n_null_fit = sc$n_null_fit, row.names = NULL)
  }))
  utils::write.csv(res, out, row.names = FALSE)
  cat(sprintf("%s rep %d done in %.1f min\n", co, r, (proc.time()[["elapsed"]] - t0) / 60))
}
cat("REPEATED_CV_DONE\n")
