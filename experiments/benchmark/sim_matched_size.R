# survfcbf top-K (q = 2) against the rankers and mRMR at a matched panel size, on the
# simulation grid with known truth. K = 3 (the number of driver clusters) and K = 10. The
# 300 datasets are those of the cesga_core grid (10 configurations, A_cens50 dropped as a
# copy of A_n600, seed 1000 + replicate), with the same five-fold harness and ridge Cox
# learner. Per fold, the selection is scored against the cluster-level truth with
# score_selection_clusterwise(): sensitivity (driver clusters recovered), FDR (distractor
# and noise selections), redundant copies of a driver cluster (parsimony_penalty) and
# exact recovery, plus the mean absolute correlation within the panel.
#
# One CSV per cell in SIM_OUTDIR (default results/sim_matched_size/cells), with the cells
# strided over N_TASKS array tasks. The published run (CESGA, 2026-11-19, about 20 min):
#   sbatch -p short --qos=short --time=03:00:00 --mem=8G --cpus-per-task=1 --array=1-60 \
#     --export=ALL,N_TASKS=60 cesga/submit_script.sbatch sim_matched_size.R
# Then, with every cell present, SIM_COMBINE=1 binds them into
# results/raw_sim_matched_size.csv, the file paper/make_statistics.R reads. MIM and mRMR
# at q = 2, the binning of top-K, were appended in a second run, whose other rows
# reproduce the first run exactly.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners); library(stabm)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
stopifnot("this needs the tie-break build of survfcbf" = "tiebreak" %in% names(formals(survfcbf_select)))
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
OUT <- Sys.getenv("SIM_OUTDIR", "results/sim_matched_size/cells")
COMBINED <- "results/raw_sim_matched_size.csv"
source("grid_lib.R")
source("realdata/uno_truncated.R")
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]

cfg <- grid_config("cesga_core")
scn <- build_scenarios(cfg)
scn <- scn[vapply(scn, function(s) s$name, "") != "A_cens50"]      # same datasets as A_n600
cells <- expand.grid(s = seq_along(scn), rep = seq_len(cfg$R))
cell_file <- function(i) file.path(OUT, sprintf("%s_rep%02d.csv", scn[[cells$s[i]]]$name, cells$rep[i]))

if (nzchar(Sys.getenv("SIM_COMBINE", ""))) {
  fs <- vapply(seq_len(nrow(cells)), cell_file, "")
  stopifnot("cells missing, run the array first" = all(file.exists(fs)))
  utils::write.csv(do.call(rbind, lapply(fs, utils::read.csv)), COMBINED, row.names = FALSE)
  cat("wrote", COMBINED, "\n"); quit(save = "no")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_tasks <- as.integer(Sys.getenv("N_TASKS", "1"))
mine <- seq.int(id, nrow(cells), by = n_tasks)

Ks <- c(3L, 10L)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")

run_cell <- function(sc, rep) {
  seed <- as.integer(1000L + rep)
  sim <- simulate_dataset(utils::modifyList(sc$params, list(seed = seed)))
  X <- sim$X; colnames(X) <- paste0("V", seq_len(ncol(X)))
  task <- surv_task(sim$time, sim$status, X)
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
  sels <- list()
  for (K in Ks) local({ k <- K
    sels[[sprintf("survfcbf_q2_K%d", k)]] <<- function(ti, st, Xs) utils::head(memo("fcbf", ti, Xs,
      function() survfcbf_select(ti, st, Xs, q = 2L, top = max(Ks))$selected), k)
    sels[[sprintf("mrmr_K%d", k)]] <<- function(ti, st, Xs) utils::head(memo("mrmr", ti, Xs,
      function() run_baseline(ent("mrmr_surv"), ti, st, Xs, prop = max(Ks) / ncol(Xs), seed = 1L)$selected), k)
    for (b in c("variance", "carss", "mim_surv")) local({ bb <- b
      sels[[sprintf("%s_K%d", bb, k)]] <<- function(ti, st, Xs) cut_ranking(rank_of(bb, ti, st, Xs), k = k) })
  })
  ## MIM and mRMR at the binning of top-K (q = 2), so the discretisation is not confounded
  ## with the pruning (by default they take q = 10)
  for (K in Ks) local({ k <- K
    sels[[sprintf("mim_surv_q2_K%d", k)]] <<- function(ti, st, Xs) cut_ranking(memo("mim_q2", ti, Xs,
      function() mim_surv_select(ti, st, Xs, q = 2L)$score), k = k)
    sels[[sprintf("mrmr_q2_K%d", k)]] <<- function(ti, st, Xs) utils::head(memo("mrmr_q2", ti, Xs,
      function() mrmr_surv_select(ti, st, Xs, prop = max(Ks) / ncol(Xs), q = 2L)$selected), k)
  })
  sim_mat <- NULL
  do.call(rbind, lapply(names(sels), function(nm) {
    bm <- nested_cv_select(task, sels[[nm]], L$clone(deep = TRUE), outer_folds = cfg$outer, seed = seed)
    scb <- score_benchmark(bm)
    if (is.null(sim_mat)) sim_mat <<- similarity_matrix(bm, scope = "full")
    st <- benchmark_stability(bm, theta = 0.9, seed = seed, sim_mat = sim_mat)
    per_fold <- lapply(bm$folds, function(fd) {
      idx <- match(as.character(fd$features), colnames(X))
      if (!length(idx)) return(NULL)
      s <- score_selection_clusterwise(idx, sim$manifest)
      r <- if (length(idx) >= 2L) { m <- abs(stats::cor(X[fd$train_set, idx, drop = FALSE])); mean(m[upper.tri(m)]) } else NA_real_
      c(sensitivity = s$sensitivity, fdr = s$fdr, redundant_copies = s$parsimony_penalty,
        exact = as.numeric(s$exact_recovery), mean_abs_r = r)
    })
    pf <- do.call(rbind, per_fold)
    data.frame(scenario = sc$name, replicate = rep, events = sum(sim$status), selector = nm,
               uno_c = scb$summary$mean[scb$summary$metric == "uno_c"],
               uno_c_trunc = mean(uno_by_fold(bm)$uno_c_trunc),
               stability = st$stability, parsimony = st$parsimony_mean,
               sensitivity = mean(pf[, "sensitivity"]), fdr = mean(pf[, "fdr"]),
               redundant_copies = mean(pf[, "redundant_copies"]), exact = mean(pf[, "exact"]),
               mean_abs_r = mean(pf[, "mean_abs_r"], na.rm = TRUE),
               n_failed = scb$n_failed, n_null_fit = scb$n_null_fit, row.names = NULL)
  }))
}

for (i in mine) {
  sc <- scn[[cells$s[i]]]; rep <- cells$rep[i]
  out <- cell_file(i)
  if (file.exists(out)) { cat(basename(out), "exists, skip\n"); next }
  t0 <- proc.time()[["elapsed"]]
  res <- run_cell(sc, rep)
  utils::write.csv(res, out, row.names = FALSE)
  cat(sprintf("%s rep %d done: %d rows, %.1f min\n", sc$name, rep, nrow(res),
              (proc.time()[["elapsed"]] - t0) / 60))
}
if (id == 1L) writeLines(utils::capture.output(utils::sessionInfo()), file.path(OUT, "sessionInfo_sim.txt"))
cat("SIM_MATCHED_SIZE_DONE\n")
