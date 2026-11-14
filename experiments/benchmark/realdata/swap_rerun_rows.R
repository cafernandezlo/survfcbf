# Swap the rows of a partial rerun into the full result tables (used for the Spearman
# tie-break of survfcbf, November 2026). A selector's row depends only on the data, the
# seed and its own code: the folds, the per-fold seeds and the stability similarity
# matrix do not depend on which other selectors run, so rerunning one selector and
# replacing its rows is equivalent to rerunning everything.
#
# The rerun includes a control, survfcbf_topk_index (top-K with the old column-order
# tie-break), which must reproduce the stored survfcbf_topk rows exactly in every cohort
# and in METABRIC. That checks the equivalence on the real data before anything is
# replaced. Runtime columns are excluded from the check, since they are wall-clock times.
#
# Inputs (from the CESGA reruns): results/realdata/cells_tiebreak/, results/
# realdata_cardinality_tiebreak/, results/realdata_metabric_tiebreak/. The tables it
# rewrites keep their previous version beside them with a _before_tiebreak suffix.
#
# Usage (from experiments/benchmark):  Rscript --no-init-file realdata/swap_rerun_rows.R
bench <- Sys.getenv("BENCH_DIR", unset = ".")
setwd(bench)
read_dir <- function(d) do.call(rbind, lapply(list.files(d, "\\.csv$", full.names = TRUE),
                                              utils::read.csv, check.names = FALSE))
keep_copy <- function(f) {
  b <- sub("\\.csv$", "_before_tiebreak.csv", f)
  if (!file.exists(b)) file.copy(f, b)
}
## the control must equal the stored rows on every non-timing column
check_control <- function(old, new, key) {
  ref <- old[old$selector == "survfcbf_topk", ]
  ctl <- new[new$selector == "survfcbf_topk_index", ]
  ref <- ref[order(ref[[key]]), ]; ctl <- ctl[order(ctl[[key]]), ]
  stopifnot("control and stored rows cover different cohorts" = identical(ref[[key]], ctl[[key]]))
  cols <- setdiff(intersect(names(ref), names(ctl)), c(key, "selector", "runtime", "runtime_se"))
  for (cc in cols) {
    a <- ref[[cc]]; b <- ctl[[cc]]
    ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < 1e-9)
    if (!all(ok)) stop(sprintf("control differs from the stored survfcbf_topk rows in %s (%s)", cc,
                               paste(ref[[key]][!ok], collapse = ", ")), call. = FALSE)
  }
  cat(sprintf("control reproduces the stored survfcbf_topk rows exactly (%d rows, %d columns)\n",
              nrow(ref), length(cols)))
}
swap <- function(old, new, key, selectors) {
  for (s in selectors) for (k in unique(new[[key]][new$selector == s])) {
    i <- which(old$selector == s & old[[key]] == k); j <- which(new$selector == s & new[[key]] == k)
    stopifnot("exactly one row to replace" = length(i) == 1L && length(j) == 1L)
    old[i, names(old)] <- new[j, names(old)]
  }
  old
}

## ---- main TCGA benchmark ----------------------------------------------------------
f <- "results/benchmark_realdata.csv"
old <- utils::read.csv(f, check.names = FALSE)
new <- read_dir("results/realdata/cells_tiebreak")
stopifnot("the rerun must cover all 11 cohorts" = length(unique(new$dataset)) == 11L)
check_control(old, new, "dataset")
keep_copy(f)
utils::write.csv(swap(old, new, "dataset", c("survfcbf_topk", "survfcbf")), f, row.names = FALSE)
cat("swapped survfcbf_topk and survfcbf rows into", f, "\n")

## ---- METABRIC ---------------------------------------------------------------------------
f <- "results/realdata_metabric/benchmark_metabric.csv"
old <- utils::read.csv(f, check.names = FALSE)
new <- utils::read.csv("results/realdata_metabric_tiebreak/benchmark_metabric.csv", check.names = FALSE)
check_control(old, new, "dataset")
keep_copy(f)
utils::write.csv(swap(old, new, "dataset", "survfcbf_topk"), f, row.names = FALSE)
cat("swapped the survfcbf_topk row into", f, "\n")

## ---- matched panel size (survfcbf rows only) ---------------------------------------------
for (nf in list.files("results/realdata_cardinality_tiebreak", "^[A-Z]+\\.csv$", full.names = TRUE)) {
  f <- file.path("results/realdata_cardinality", basename(nf))
  old <- utils::read.csv(f, check.names = FALSE); new <- utils::read.csv(nf, check.names = FALSE)
  stopifnot("the cardinality rerun holds survfcbf rows only" = all(new$base == "survfcbf"),
            "different cohort" = identical(unique(old$cohort), unique(new$cohort)))
  old$key <- paste(old$base, old$K); new$key <- paste(new$base, new$K)
  for (k in new$key) old[old$key == k, setdiff(names(old), "key")] <- new[new$key == k, setdiff(names(old), "key")]
  old$key <- NULL
  b <- file.path("results/realdata_cardinality_before_tiebreak", basename(f))
  dir.create(dirname(b), showWarnings = FALSE)
  if (!file.exists(b)) file.copy(f, b)
  utils::write.csv(old, f, row.names = FALSE)
}
cat("swapped survfcbf rows into results/realdata_cardinality/*.csv\n")
cat("SWAP_DONE\n")
