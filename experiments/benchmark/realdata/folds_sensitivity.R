# Is the between-cohort spread irreducible, or would better evaluation shrink it?
#
# The required-datasets figures rest on the standard deviation of the paired per-cohort
# differences. That spread mixes genuine between-cohort heterogeneity, which only more
# datasets can address, with residual noise from estimating each cohort's performance,
# which more folds would reduce. The stored summaries cannot separate the two because the
# harness keeps a mean and a Monte Carlo SE per cell rather than per-fold scores.
#
# This answers the question empirically instead: re-run the cardinality comparison at the
# decisive K with twice the outer folds and compare the resulting paired SD against the
# five-fold run. If doubling the folds barely moves it, the spread is dominated by
# heterogeneity and no evaluation effort removes the requirement. If it drops materially,
# the required sample sizes are an overestimate and must be qualified.
#
# SLURM array, one cohort per task, coding-only, K fixed at the value where the paired
# comparison was decisive.
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
FOLDS <- as.integer(Sys.getenv("RD_FOLDS", "10"))     # the five-fold run already exists
K <- 20L

outdir <- "results/realdata_folds"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, sprintf("%s_f%d.csv", co, FOLDS))
if (file.exists(out)) { cat(co, "exists, skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
cat(sprintf("%s: n=%d ev=%d p=%d, %d outer folds, K=%d\n",
            co, length(d$time), sum(d$status), ncol(X), FOLDS, K))

cache <- new.env(parent = emptyenv())
memo <- function(key, ti, X, fn) {
  k <- sprintf("%s_%d_%.4f", key, nrow(X), sum(ti))
  if (!is.null(cache[[k]])) return(cache[[k]])
  v <- fn(); assign(k, v, envir = cache); v
}
sels <- list(
  survfcbf = function(ti, st, Xs) memo("fcbf", ti, Xs, function()
    survfcbf_select(ti, st, Xs, q = 2L, top = K)$selected),
  variance = function(ti, st, Xs) cut_ranking(apply(Xs, 2L, stats::var), k = K),
  carss    = function(ti, st, Xs) cut_ranking(
    abs(run_baseline(ent("carss"), ti, st, Xs, prop = NULL, seed = 1L)$score), k = K),
  mim_surv = function(ti, st, Xs) cut_ranking(
    run_baseline(ent("mim_surv"), ti, st, Xs, prop = NULL, seed = 1L)$score, k = K))

task <- surv_task(d$time, d$status, X)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")   # 1se -> null model (?score_benchmark)
tab <- benchmark_selectors(task, sels, L$clone(deep = TRUE), outer_folds = FOLDS,
                           manifest = NULL, seed = 1L)
res <- data.frame(cohort = co, folds = FOLDS, K = K, n = length(d$time),
                  events = sum(d$status), p = ncol(X), tab, row.names = NULL,
                  check.names = FALSE)
utils::write.csv(res, out, row.names = FALSE)
cat(sprintf("%s done (%d folds)\n", co, FOLDS))
cat("FOLDS_DONE\n")
