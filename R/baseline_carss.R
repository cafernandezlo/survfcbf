#' carss IPCW-weighted correlation selection (baseline)
#'
#' Wraps the CAR survival score of Bommert's near-best filter:
#' `carSurv::carSurvScore` computes an IPCW-weighted, decorrelated
#' (Correlation-Adjusted Regression) score per feature, and
#' `carSurv::carVarSelect` turns those scores into a self-selected set by a
#' false-non-discovery-rate rule. This is the key predictive comparator, since the
#' carss family was near the top on the integrated Brier score in Bommert et al.
#' (2022).
#'
#' @details
#' The CAR scores are signed, so the per-feature importance returned in `score` is
#' their magnitude (the natural ranking signal), while the signed scores are kept
#' in `raw_score` and are what `carVarSelect` consumes. IPCW weights are capped at
#' `maxIPCweight` (carSurv's own default of 10). The score and the censoring
#' weights are computed on the supplied (training) fold only, so the selection
#' carries no leakage. With `prop` given the top-`prop` features by importance are
#' returned, otherwise carss self-selects via `carVarSelect`. Note that
#' `carVarSelect`'s `fndr` rule (its only non-erroring method in carSurv 1.0)
#' returns the empty set at its default threshold on weak-to-moderate signal, so the
#' benchmark registry wires carss as a ranker at `prop`, faithful to its use in
#' Bommert et al. (2022), rather than via this self-selection. Deterministic given
#' the data.
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. `NULL`
#'   (default) uses carss's own `carVarSelect` set.
#' @param maxIPCweight Cap on the inverse-probability-of-censoring weights, passed
#'   to `carSurv::carSurvScore`. Default `10`.
#' @param threshold False-non-discovery-rate threshold for `carVarSelect` when
#'   `prop` is `NULL`. Default `0.05`.
#' @param seed Unused (the score is deterministic); present for interface
#'   uniformity.
#'
#' @return A list with `score` (per-feature importance, the absolute CAR score,
#'   length `p`), `selected` (the `carVarSelect` set, or the top-`prop` indices),
#'   and `raw_score` (the signed CAR scores).
#'
#' @seealso [variance_select()], [cut_ranking()]
#' @export
carss_select <- function(time, status, X, prop = NULL,
                         maxIPCweight = 10, threshold = 0.05, seed = NULL) {
  if (!requireNamespace("carSurv", quietly = TRUE)) {
    stop("carss_select requires the 'carSurv' package; install it with install.packages('carSurv').",
         call. = FALSE)
  }
  stopifnot(
    "`X` must be a numeric matrix"                 = is.matrix(X) && is.numeric(X),
    "`X`, `time`, `status` must share nrow/length" = nrow(X) == length(time) && length(time) == length(status)
  )
  if (is.logical(status)) status <- as.integer(status)

  raw <- carSurv::carSurvScore(obsTime = time, obsEvent = status, X = X, maxIPCweight = maxIPCweight)
  score <- abs(raw)                                  # CAR scores are signed; importance is magnitude
  selected <- if (is.null(prop)) {
    as.integer(carSurv::carVarSelect(raw, method = "fndr", threshold = threshold))
  } else {
    cut_ranking(score, prop = prop)
  }
  list(score = score, selected = selected, raw_score = raw)
}
