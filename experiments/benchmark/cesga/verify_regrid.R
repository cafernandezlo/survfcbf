# Verify the corrected-recovery re-run of the simulation grid against the original.
#
# The recovery columns of the June 2026 cesga_core run were computed by a
# benchmark_recovery() that scored a fold's mlr3 (alphabetically ordered) feature
# indices as manifest row numbers. The fix resolves feature NAMES instead. Because
# the fix touches only the scoring, and because selection is deterministic given the
# cell seed and independent of the downstream learner, a re-run must reproduce the
# original selections EXACTLY. That gives a strong self-check:
#
#   uno_c, harrell_c, ibs, stability, parsimony, n_failed, events  -> must be identical
#   sensitivity, fdr, exact_recovery                               -> must change
#   runtime                                                        -> free to differ (node speed)
#
# If an invariant column moves, the determinism argument is broken and the recovery
# numbers must NOT be merged into the old table: something other than the scorer changed.
#
# Usage (from experiments/benchmark):
#   Rscript --no-init-file cesga/verify_regrid.R <old_raw.csv> <new_raw.csv>

args <- commandArgs(trailingOnly = TRUE)
old_f <- if (length(args) >= 1L) args[1] else
  "results/raw_cesga_core_buggy_recovery_jun2026.csv"
new_f <- if (length(args) >= 2L) args[2] else "results/raw_cesga_core.csv"

stopifnot("old raw file not found" = file.exists(old_f),
          "new raw file not found" = file.exists(new_f))
old <- utils::read.csv(old_f, stringsAsFactors = FALSE)
new <- utils::read.csv(new_f, stringsAsFactors = FALSE)

KEYS      <- c("scenario", "map", "learner", "replicate", "selector")
INVARIANT <- c("events", "uno_c", "harrell_c", "ibs", "stability", "parsimony", "n_failed")
CHANGED   <- c("sensitivity", "fdr", "exact_recovery")

cat(sprintf("old: %s  (%d rows)\nnew: %s  (%d rows)\n\n",
            old_f, nrow(old), new_f, nrow(new)))

## ---- coverage ------------------------------------------------------------------
kold <- do.call(paste, c(old[KEYS], sep = "|"))
knew <- do.call(paste, c(new[KEYS], sep = "|"))
miss <- setdiff(kold, knew)
extra <- setdiff(knew, kold)
cat(sprintf("cells in old missing from new: %d\ncells in new missing from old: %d\n\n",
            length(miss), length(extra)))
if (length(miss)) cat("  missing (first 10):\n    ",
                      paste(utils::head(miss, 10), collapse = "\n    "), "\n\n")

m <- merge(old, new, by = KEYS, suffixes = c(".old", ".new"))
cat(sprintf("matched rows: %d\n\n", nrow(m)))
if (!nrow(m)) stop("no matched rows: check the keys")

## ---- the columns that must NOT move ---------------------------------------------
cat("=== INVARIANTS (must be identical) ===\n")
fails <- character(0)
for (v in INVARIANT) {
  a <- m[[paste0(v, ".old")]]; b <- m[[paste0(v, ".new")]]
  if (is.null(a) || is.null(b)) { cat(sprintf("  %-12s MISSING from one of the two tables\n", v)); next }
  both_na <- is.na(a) & is.na(b)
  one_na  <- xor(is.na(a), is.na(b))
  d <- abs(a - b); d[both_na] <- 0
  worst <- suppressWarnings(max(d, na.rm = TRUE))
  if (!is.finite(worst)) worst <- 0
  bad <- sum(one_na, na.rm = TRUE) + sum(d > 1e-8, na.rm = TRUE)
  status <- if (bad == 0L) "OK" else "DIFFERS"
  if (bad > 0L) fails <- c(fails, v)
  cat(sprintf("  %-12s %-8s  max|diff| = %.3e   differing rows = %d  (unpaired NA = %d)\n",
              v, status, worst, bad, sum(one_na, na.rm = TRUE)))
}

## ---- the columns that SHOULD move ------------------------------------------------
cat("\n=== RECOVERY (must change) ===\n")
for (v in CHANGED) {
  a <- m[[paste0(v, ".old")]]; b <- m[[paste0(v, ".new")]]
  if (is.null(a) || is.null(b)) { cat(sprintf("  %-15s MISSING\n", v)); next }
  ok <- is.finite(a) & is.finite(b)
  moved <- sum(ok & abs(a - b) > 1e-8)
  cat(sprintf("  %-15s old %.4f -> new %.4f   (mean over %d valid rows, %d rows change, NA: old %d / new %d)\n",
              v, mean(a[ok]), mean(b[ok]), sum(ok), moved, sum(is.na(a)), sum(is.na(b))))
}

## ---- per-selector recovery shift --------------------------------------------------
cat("\n=== sensitivity by selector (mean over valid replicates) ===\n")
for (s in sort(unique(m$selector))) {
  z <- m[m$selector == s, ]
  ok <- is.finite(z$sensitivity.old) & is.finite(z$sensitivity.new)
  if (!any(ok)) { cat(sprintf("  %-16s no valid rows\n", s)); next }
  d <- z$sensitivity.new[ok] - z$sensitivity.old[ok]
  cat(sprintf("  %-16s old %.3f -> new %.3f   diff %+.3f +/- %.3f (n=%d)\n",
              s, mean(z$sensitivity.old[ok]), mean(z$sensitivity.new[ok]),
              mean(d), stats::sd(d) / sqrt(length(d)), length(d)))
}

## ---- verdict -----------------------------------------------------------------------
cat("\n=== VERDICT ===\n")
if (length(fails)) {
  cat("FAIL: these invariants moved ->", paste(fails, collapse = ", "), "\n")
  cat("The rerun did NOT reproduce the original selections. Do not merge: find what else changed\n")
  cat("(package version in the container, seed, feature order, tier config).\n")
  quit(status = 1L)
}
cat("OK: every invariant reproduces, so the rerun repeated the original selections exactly\n")
cat("and the only difference is the corrected scorer.\n")
