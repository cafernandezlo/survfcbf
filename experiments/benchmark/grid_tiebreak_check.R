# Does the Spearman tie-break change the survfcbf selections of the simulation grid? The
# grid ran before the tie-break existed, with ties in relevance broken by column order.
# For 3 replicates of each of the 10 cesga_core configurations (30 datasets, A_cens50
# dropped as a copy of A_n600, seed 1000 + replicate) and each of the 5 outer training
# sets, the permutation mode is run with tiebreak = "spearman" and then with
# tiebreak = "index" (the column-order rule) at the SAME delta, for the settings the grid
# used: the default q (Bommert rule, B = 200) and q = 4, 7, 10 (B = 100, the tuned-q
# grid). Folds and per-fold seeds as in the benchmark harness.
#
# Writes results/grid_tiebreak_check.csv, read by paper/make_supplementary_tables.R. The
# published run (local, 2026-11-20, 14 workers, 75 min) loaded the package from the source
# tree with SURVFCBF_SRC=../.. and gave the same set in the same order in all 600
# comparisons, with no tie in relevance above delta.
SRC <- Sys.getenv("SURVFCBF_SRC", "")
load_pkg <- function() {
  if (nzchar(SRC)) suppressMessages(devtools::load_all(SRC, quiet = TRUE))
  else suppressMessages(library(survfcbf))
}
load_pkg()
stopifnot("this needs the tie-break build of survfcbf" = "tiebreak" %in% names(formals(survfcbf_select)))
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
OUT <- "results/grid_tiebreak_check.csv"
WORKERS <- as.integer(Sys.getenv("N_WORKERS", "14"))
source("grid_lib.R")
cfg <- grid_config("cesga_core")
scn <- build_scenarios(cfg)
scn <- scn[vapply(scn, function(s) s$name, "") != "A_cens50"]
units <- expand.grid(s = seq_along(scn), rep = 1:3)
SET <- list(default = list(q = NULL, B = 200L), q4 = list(q = 4L, B = 100L),
            q7 = list(q = 7L, B = 100L), q10 = list(q = 10L, B = 100L))
one <- function(i) {
  load_pkg()
  sc <- scn[[units$s[i]]]; rep <- units$rep[i]; seed <- as.integer(1000L + rep)
  sim <- simulate_dataset(utils::modifyList(sc$params, list(seed = seed)))
  folds <- survfcbf:::stratified_folds(sim$status, cfg$outer, seed)
  out <- list()
  for (f in seq_len(cfg$outer)) {
    tr <- folds != f
    fs <- as.integer((seed * 1000 + f) %% 2147483647)
    for (nm in names(SET)) {
      st <- SET[[nm]]
      a <- survfcbf:::with_seed(fs, suppressWarnings(survfcbf_select(sim$time[tr], sim$status[tr], sim$X[tr, , drop = FALSE],
                                                 q = st$q, B = st$B, tiebreak = "spearman")))
      b <- survfcbf:::with_seed(fs, suppressWarnings(survfcbf_select(sim$time[tr], sim$status[tr], sim$X[tr, , drop = FALSE],
                                                 q = st$q, B = st$B, delta = a$delta, tiebreak = "index")))
      rel <- a$relevance; above <- rel[rel >= a$delta]
      out[[length(out) + 1L]] <- data.frame(
        scenario = sc$name, replicate = rep, fold = f, setting = nm, delta = a$delta,
        n_relevant = length(above), n_tied_relevant = sum(duplicated(above) | duplicated(above, fromLast = TRUE)),
        n_sel_spearman = length(a$selected), n_sel_index = length(b$selected),
        same_set = setequal(a$selected, b$selected),
        same_order = identical(as.integer(a$selected), as.integer(b$selected)))
    }
  }
  do.call(rbind, out)
}
cl <- parallel::makeCluster(WORKERS)
parallel::clusterExport(cl, c("SRC", "load_pkg", "cfg", "scn", "units", "SET"))
t0 <- proc.time()[["elapsed"]]
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(units)), one))
parallel::stopCluster(cl)
utils::write.csv(res, OUT, row.names = FALSE)
cat(sprintf("done in %.1f min, %d comparisons\n", (proc.time()[["elapsed"]] - t0) / 60, nrow(res)))
print(aggregate(cbind(same_set, same_order, n_tied_relevant > 0, n_sel_spearman != n_sel_index) ~ setting, res, sum))
