# SLURM array driver for the Bommert real-data benchmark. Task i
# loads dataset i (load_bommert, from the pre-downloaded ARFF cache -- compute nodes
# have no internet, so the loader resolves cached files offline by id), runs the
# nested-CV harness, and writes results/realdata/cells/<name>.csv. There is no known
# truth on real data, so manifest = NULL: prediction (Uno C / IBS), stability and
# parsimony only -- recovery is assessed later by biological enrichment, not here.
# Run inside the bioconductor R 4.5.2 apptainer container (see submit_realdata.sbatch).
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners)
  library(stabm); library(jsonlite); library(curl); library(farff)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = "")            # apptainer cwd is unreliable
if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R")

# reduced rigor for the p~23k regime; q_grid includes small q (the events lever).
# survfcbf_topk is the explicit-cardinality mode (q=2 top-K): the redesign
# fix for the delta self-selection collapse, the one to validate against survfcbf here.
cfg <- list(outer = 5L, inner_folds = 3L, B = 100L, q_grid = c(2L, 4L, 7L, 10L), prop = 0.1, topk = 10L)
KEEP <- c("survfcbf", "survfcbf_topk", "survfcbf_cvq", "variance", "carss", "cox_glmnet_l1", "mim_surv")
if (nzchar(Sys.getenv("RD_SMOKE", ""))) {           # fast end-to-end smoke: reduced rigor, 4 selectors
  cfg <- list(outer = 2L, inner_folds = 2L, B = 20L, q_grid = c(4L, 7L), prop = 0.1, topk = 10L)
  KEEP <- c("variance", "survfcbf_topk", "survfcbf_cvq", "carss")
}
## Partial reruns. RD_KEEP (comma-separated selector names) runs a subset and RD_CELLS
## writes its cells to another directory, so a rerun never overwrites the full run. A
## selector's row does not depend on which others run with it: the folds, the per-fold
## seeds and the stability similarity matrix depend only on the data and the seed.
## "survfcbf_topk_index" is top-K with the old column-order tie-break, a control that
## must reproduce the stored survfcbf_topk rows exactly.
if (nzchar(Sys.getenv("RD_KEEP", ""))) KEEP <- strsplit(Sys.getenv("RD_KEEP"), ",", fixed = TRUE)[[1]]
cells_dir <- Sys.getenv("RD_CELLS", "results/realdata/cells")

datasets <- names(BOMMERT_IDS)                           # 11 cohorts
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot("SLURM_ARRAY_TASK_ID out of range" = !is.na(id) && id >= 1L && id <= length(datasets))
ds <- datasets[id]

t0 <- proc.time()[["elapsed"]]
d <- load_bommert(ds)                                   # (time, status, X) from cache
task <- surv_task(d$time, d$status, d$X)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")   # 1se -> null model (?score_benchmark)
sels <- default_selectors(prop = cfg$prop, learner = L, q_grid = cfg$q_grid,
                          cv_args = list(inner_folds = cfg$inner_folds, B = cfg$B), topk = cfg$topk)
sels[["survfcbf_topk_index"]] <- function(ti, st, X)
  survfcbf_select(ti, st, X, q = 2L, top = cfg$topk, tiebreak = "index")$selected
sels <- sels[intersect(KEEP, names(sels))]
stopifnot("RD_KEEP names no known selector" = length(sels) > 0L)
tab <- benchmark_selectors(task, sels, L$clone(deep = TRUE), outer_folds = cfg$outer,
                           manifest = NULL, seed = 1L)
res <- data.frame(dataset = ds, n = length(d$time), events = sum(d$status), p = ncol(d$X),
                  tab, row.names = NULL, check.names = FALSE)
outdir <- cells_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(outdir, sprintf("%s.csv", ds)), row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(outdir, sprintf("sessionInfo_%s.txt", ds)))
cat(sprintf("%s done: n=%d events=%d p=%d -> %d rows, %.1f min\n",
            ds, length(d$time), sum(d$status), ncol(d$X), nrow(res), (proc.time()[["elapsed"]] - t0) / 60))
