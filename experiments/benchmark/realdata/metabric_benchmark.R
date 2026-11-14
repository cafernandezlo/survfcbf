# METABRIC external-cohort benchmark (external validation). Runs the
# same leakage-safe nested-CV harness on the METABRIC breast cohort (n~1979, events~1143,
# p~19.8k complete profile, log2 microarray, OS) with survfcbf_topk vs the ranking baselines. survfcbf_cvq
# is dropped here (the expensive continuity selector); the external validation is of the
# redesign (top-K) against variance/carss/mim/cox-lasso. METABRIC has ~1143 events (vs the
# 52-152 of the Bommert cohorts), so this is also the high-event corner of the regime.
# METABRIC_DIR points at the .rds data (set in the sbatch). Run in the apptainer container.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners); library(stabm)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench)) setwd(bench)
source("realdata/load_metabric.R")

cfg  <- list(outer = 5L, inner = 3L, B = 100L, q_grid = c(2L, 4L, 7L, 10L), prop = 0.1, topk = 10L)
KEEP <- c("survfcbf_topk", "variance", "carss", "cox_glmnet_l1", "mim_surv")
## partial rerun: RD_KEEP runs a subset, RD_OUTFILE writes it elsewhere (see run_dataset.R)
if (nzchar(Sys.getenv("RD_KEEP", ""))) KEEP <- strsplit(Sys.getenv("RD_KEEP"), ",", fixed = TRUE)[[1]]
OUTFILE <- Sys.getenv("RD_OUTFILE", "results/realdata_metabric/benchmark_metabric.csv")

t0 <- proc.time()[["elapsed"]]
d  <- load_metabric()
cat(sprintf("METABRIC: n=%d events=%d p=%d\n", length(d$time), sum(d$status), ncol(d$X)))
task <- surv_task(d$time, d$status, d$X)
L    <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")   # 1se -> null model (?score_benchmark)
sels <- default_selectors(prop = cfg$prop, topk = cfg$topk)     # no learner -> no cvq/atprop
sels[["survfcbf_topk_index"]] <- function(ti, st, X)           # control: old column-order tie-break
  survfcbf_select(ti, st, X, q = 2L, top = cfg$topk, tiebreak = "index")$selected
sels <- sels[intersect(KEEP, names(sels))]
tab  <- benchmark_selectors(task, sels, L$clone(deep = TRUE), outer_folds = cfg$outer,
                            manifest = NULL, seed = 1L)
res  <- data.frame(dataset = "METABRIC", n = length(d$time), events = sum(d$status),
                   p = ncol(d$X), tab, row.names = NULL, check.names = FALSE)
dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, OUTFILE, row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), sub("\\.csv$", "_sessionInfo.txt", OUTFILE))
cat(sprintf("METABRIC benchmark done: %d rows, %.1f min\n", nrow(res), (proc.time()[["elapsed"]] - t0) / 60))
cat("METABRIC_BENCH_DONE\n")
