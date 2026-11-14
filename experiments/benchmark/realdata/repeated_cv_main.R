# Repeated cross-validation of the headline comparison of the Application Note: top-K
# (q = 2, ten genes) against the three rankers at the conventional ten-percent cut-off, on
# the full feature set of the main benchmark (RNA and miRNA, as run_dataset.R). The
# published benchmark uses one five-fold split per cohort (seed 1). This repeats the split
# ten times (seeds 1 to 10, seed 1 being the published split, which it must reproduce) and
# also scores Uno's C truncated at the 80th percentile of each training fold's follow-up
# (realdata/uno_truncated.R).
#
# SLURM array over (cohort, repetition) units, with the repetitions of a cohort kept
# together so the data and the similarity matrix (about 4 GB at p ~ 23k) load once per task:
#   sbatch -p short --qos=short --time=06:00:00 --mem=64G --cpus-per-task=2 --array=1-22 \
#     --export=ALL,N_TASKS=22 cesga/submit_script.sbatch realdata/repeated_cv_main.R
# Writes <cohort>_rep<r>.csv to RD_OUTDIR (default results/realdata_repeated_cv_main).
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
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_repeated_cv_main")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

REPS <- as.integer(Sys.getenv("RCV_REPS", "10"))
units <- expand.grid(rep = seq_len(REPS), cohort = names(BOMMERT_IDS), stringsAsFactors = FALSE)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_tasks <- as.integer(Sys.getenv("N_TASKS", as.character(nrow(units))))
per <- ceiling(nrow(units) / n_tasks)
mine <- seq.int((id - 1L) * per + 1L, min(id * per, nrow(units)))
PROP <- 0.1; TOPK <- 10L                                 # as cfg in run_dataset.R
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")
## the selectors exactly as default_selectors() builds them for run_dataset.R
sels <- list(
  survfcbf_topk = function(ti, st, X) survfcbf_select(ti, st, X, q = 2L, top = TOPK)$selected,
  variance = function(ti, st, X) run_baseline(ent("variance"), ti, st, X, prop = PROP)$selected,
  carss    = function(ti, st, X) run_baseline(ent("carss"), ti, st, X, prop = PROP)$selected,
  mim_surv = function(ti, st, X) run_baseline(ent("mim_surv"), ti, st, X, prop = PROP)$selected)

loaded <- ""
for (i in mine) {
  co <- units$cohort[i]; r <- units$rep[i]
  out <- file.path(outdir, sprintf("%s_rep%02d.csv", co, r))
  if (file.exists(out)) { cat(basename(out), "exists, skip\n"); next }
  t0 <- proc.time()[["elapsed"]]
  if (loaded != co) {
    d <- load_bommert(co)
    task <- surv_task(d$time, d$status, d$X)
    sim_mat <- NULL                         # outcome-free and fold-free: one per cohort
    loaded <- co
  }
  res <- do.call(rbind, lapply(names(sels), function(nm) {
    bm <- nested_cv_select(task, sels[[nm]], L$clone(deep = TRUE), outer_folds = 5L, seed = r)
    sc <- score_benchmark(bm)
    if (is.null(sim_mat)) sim_mat <<- similarity_matrix(bm, scope = "full")
    st <- benchmark_stability(bm, theta = 0.9, seed = r, sim_mat = sim_mat)
    s <- stats::setNames(sc$summary$mean, sc$summary$metric)
    uf <- uno_by_fold(bm)
    data.frame(cohort = co, repetition = r, seed = r, p = ncol(d$X), selector = nm,
               uno_c = s[["uno_c"]], uno_c_trunc = mean(uf$uno_c_trunc),
               harrell_c = s[["harrell_c"]], ibs = s[["ibs"]],
               stability = st$stability, parsimony = st$parsimony_mean,
               n_failed = sc$n_failed, n_null_fit = sc$n_null_fit, row.names = NULL)
  }))
  utils::write.csv(res, out, row.names = FALSE)
  cat(sprintf("%s rep %d done in %.1f min\n", co, r, (proc.time()[["elapsed"]] - t0) / 60))
}
cat("REPEATED_CV_MAIN_DONE\n")
