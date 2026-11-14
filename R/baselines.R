#' Cut a per-feature ranking to a top set
#'
#' Turns a per-feature importance `score` (higher meaning more relevant) into the
#' index set of the top features, the shared bridge that lets the benchmark
#' evaluate ranking baselines at FCBF's own cardinality or at inner-CV-tuned
#' thresholds. Supply either a proportion `prop` of features to
#' keep or an explicit count `k`.
#'
#' @details
#' Ties are broken by ascending feature index, so the cut is deterministic. With
#' `prop` the kept count is `max(1, round(prop * p))`; with `k` it is `min(k, p)`.
#'
#' @param score Numeric per-feature importance, length `p` (no `NA`).
#' @param prop Proportion of features to keep, in `(0, 1]`. Give this or `k`.
#' @param k Number of features to keep, a positive integer. Give this or `prop`.
#'
#' @return Integer feature indices of the top set, ordered by decreasing score.
#'
#' @seealso [variance_select()]
#' @export
cut_ranking <- function(score, prop = NULL, k = NULL) {
  stopifnot("`score` must be a numeric vector without NA" = is.numeric(score) && length(score) >= 1L && !anyNA(score))
  p <- length(score)
  if (is.null(k)) {
    if (is.null(prop)) stop("provide either `prop` or `k`.", call. = FALSE)
    stopifnot("`prop` must be a single number in (0, 1]" =
                length(prop) == 1L && is.numeric(prop) && prop > 0 && prop <= 1)
    k <- max(1L, round(prop * p))
  } else {
    stopifnot("`k` must be a single positive integer" =
                length(k) == 1L && is.numeric(k) && is.finite(k) && k >= 1)
    k <- min(as.integer(k), p)
  }
  ## top-k by score descending, ties broken by ascending index (deterministic)
  order(score, seq_len(p), decreasing = c(TRUE, FALSE), method = "radix")[seq_len(k)]
}

#' Marginal-variance feature selection (sanity-floor baseline)
#'
#' Ranks features by their marginal variance, the trivial filter that beat the
#' mutual-information family on the integrated Brier score in Bommert et al.
#' (2022) and so is the baseline `survfcbf` must be honest about.
#' It is unsupervised: `time` and `status` are accepted for a common selector
#' interface but ignored.
#'
#' @param time,status Right-censored outcome (ignored; present for interface
#'   uniformity with the supervised baselines).
#' @param X Numeric `n` by `p` feature matrix.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. When `NULL`
#'   (default) the function is a pure ranker and `selected` is `NULL`.
#' @param seed Unused (no randomness); present for interface uniformity.
#'
#' @return A list with `score` (per-feature variance, length `p`) and `selected`
#'   (top-`prop` indices, or `NULL` when `prop` is not given).
#'
#' @seealso [cut_ranking()], [nofilter_select()]
#' @export
variance_select <- function(time, status, X, prop = NULL, seed = NULL) {
  stopifnot("`X` must be a numeric matrix" = is.matrix(X) && is.numeric(X))
  score <- apply(X, 2L, stats::var)
  selected <- if (is.null(prop)) NULL else cut_ranking(score, prop = prop)
  list(score = score, selected = selected)
}

#' No-filter pass-through (sanity-floor baseline)
#'
#' Keeps every feature, the "no selection" reference point for the benchmark.
#' The `score` is constant (it carries no ranking) and
#' `selected` is all `p` features regardless of `prop`.
#'
#' @param time,status Right-censored outcome (ignored; present for interface
#'   uniformity).
#' @param X Numeric `n` by `p` feature matrix.
#' @param prop,seed Unused; present for interface uniformity.
#'
#' @return A list with `score` (a length-`p` vector of zeros) and `selected` (all
#'   `p` feature indices).
#'
#' @seealso [variance_select()]
#' @export
nofilter_select <- function(time, status, X, prop = NULL, seed = NULL) {
  stopifnot("`X` must be a numeric matrix" = is.matrix(X) && is.numeric(X))
  p <- ncol(X)
  list(score = rep(0, p), selected = seq_len(p))
}

#' Registry of baseline selectors for the benchmark
#'
#' Lists every baseline with the metadata the benchmark needs to
#' iterate over them uniformly: the selector function, whether it
#' self-selects a set or is a ranking that must be cut at a cardinality, and the
#' Suggested package it requires.
#'
#' @details
#' `type` is `"set"` for selectors that return their own `selected` with
#' `prop = NULL` (no-filter, the lasso) and `"ranker"` for those that need a `prop`
#' to produce a set (variance, carss, the MI family, ridge, RSF). carss is ranked at
#' `prop` (faithful to Bommert et al. 2022) because `carVarSelect`'s self-selection
#' degenerates to the empty set at its default `fndr` threshold. `self_selects` is
#' the operational flag (whether `selected` is non-`NULL` at `prop = NULL`). Entries
#' that fix a configuration of a shared selector (the L1 and L2 Cox-glmnet) carry it
#' in `args`. `package` is the Suggested dependency (`NA` for base-only selectors),
#' so a caller can drop entries whose package is absent. Use [run_baseline()] to
#' invoke an entry. The survival FCBF method itself ([survfcbf_select()]) is the
#' subject of the benchmark, not a registry baseline.
#'
#' @return A list of records, each `name`, `fn` (the selector), `type`,
#'   `self_selects`, `package`, and `args` (fixed extra arguments).
#'
#' @seealso [run_baseline()], [variance_select()], [survfcbf_select()]
#' @export
baseline_registry <- function() {
  list(
    list(name = "variance",      fn = variance_select,   type = "ranker", self_selects = FALSE, package = NA_character_, args = list()),
    list(name = "nofilter",      fn = nofilter_select,   type = "set",    self_selects = TRUE,  package = NA_character_, args = list()),
    list(name = "mim_surv",      fn = mim_surv_select,   type = "ranker", self_selects = FALSE, package = "praznik",     args = list()),
    list(name = "mrmr_surv",     fn = mrmr_surv_select,  type = "ranker", self_selects = FALSE, package = "praznik",     args = list()),
    list(name = "carss",         fn = carss_select,      type = "ranker", self_selects = FALSE, package = "carSurv",     args = list()),
    list(name = "cox_glmnet_l1", fn = cox_glmnet_select, type = "set",    self_selects = TRUE,  package = "glmnet",      args = list(alpha = 1)),
    list(name = "cox_glmnet_l2", fn = cox_glmnet_select, type = "ranker", self_selects = FALSE, package = "glmnet",      args = list(alpha = 0)),
    list(name = "rsf",           fn = rsf_select,        type = "ranker", self_selects = FALSE, package = "ranger",      args = list()),
    list(name = "grplasso_ensemble", fn = grplasso_ensemble_select, type = "set", self_selects = TRUE, package = "grpreg", args = list())
  )
}

#' Run a baseline-registry entry on a training fold
#'
#' Invokes a [baseline_registry()] entry leakage-free, threading `prop` and `seed`
#' and applying the entry's fixed `args`.
#'
#' @param entry A registry record from [baseline_registry()].
#' @param time,status,X Training-fold outcome and feature matrix.
#' @param prop Optional cardinality proportion for ranking selectors.
#' @param seed Optional seed for the stochastic selectors.
#'
#' @return The selector's `list(score, selected, ...)`.
#' @seealso [baseline_registry()]
#' @export
run_baseline <- function(entry, time, status, X, prop = NULL, seed = NULL) {
  stopifnot("`entry` must be a baseline registry record with a function `fn`" =
              is.list(entry) && is.function(entry$fn))
  do.call(entry$fn, c(list(time, status, X, prop = prop, seed = seed),
                      if (is.null(entry$args)) list() else entry$args))
}
