#' Build the survival linear predictor from features and the truth manifest
#'
#' Forms the per-observation linear predictor `eta` that drives the hazard,
#' combining the additive effects of the signal (driver) features with a single
#' standardised interaction product term. `beta_scale` is the global
#' relevance-strength knob.
#'
#' @details
#' The additive part is `X[, signal] %*% (beta_scale * beta[signal])`, where the
#' signal features are those with non-zero `beta` in the manifest (the drivers).
#' If the manifest carries an interaction pair (two features with
#' `role == "interaction"`) their standardised product adds
#' `beta_scale * beta_int * std(Xa) * std(Xb)`. Both interaction features are
#' standardised so that, being independent, their product has mean zero: the
#' marginal association of either feature with `eta` is then near zero while the
#' joint (product) association is not, the limitation that FCBF's univariate
#' relevance leg is meant to exhibit.
#'
#' Both terms carry `beta_scale`, so `eta` is constant (zero) when
#' `beta_scale = 0`. The standardisation here is part of the data-generating
#' process, not estimation, so it raises no train/test leakage concern.
#'
#' @param X Numeric `n` by `p` feature matrix, e.g. `X` from [sample_features()].
#' @param manifest The data.frame from [make_sigma_blocks()], with `role` and
#'   `beta` columns and one row per column of `X`.
#' @param beta_scale Numeric scalar. Global relevance-strength multiplier applied
#'   to every effect.
#' @param beta_int Numeric scalar. Relative strength of the interaction product
#'   term, in the same units as the driver betas. Default `1`.
#'
#' @return A numeric vector `eta` of length `nrow(X)`.
#'
#' @seealso [make_sigma_blocks()], [sample_features()]
#' @export
build_linear_predictor <- function(X, manifest, beta_scale, beta_int = 1) {
  stopifnot(
    "`X` must be a numeric matrix"                   = is.matrix(X) && is.numeric(X),
    "`manifest` must have `role` and `beta` columns" = is.data.frame(manifest) && all(c("role", "beta") %in% names(manifest)),
    "`manifest` must have one row per column of `X`" = nrow(manifest) == ncol(X),
    "`beta_scale` must be a single finite number"    = length(beta_scale) == 1L && is.numeric(beta_scale) && is.finite(beta_scale),
    "`beta_int` must be a single finite number"      = length(beta_int) == 1L && is.numeric(beta_int) && is.finite(beta_int)
  )

  n <- nrow(X)
  eta <- numeric(n)

  ## additive driver (signal) effects
  signal <- which(manifest$beta != 0)
  if (length(signal) > 0L) {
    eta <- eta + as.vector(X[, signal, drop = FALSE] %*% (beta_scale * manifest$beta[signal]))
  }

  ## interaction product term (standardised pair), also scaled by beta_scale so
  ## that eta vanishes when beta_scale == 0
  inter <- which(manifest$role == "interaction")
  if (length(inter) == 2L) {
    eta <- eta + beta_scale * beta_int * std_unit(X[, inter[1L]]) * std_unit(X[, inter[2L]])
  } else if (length(inter) != 0L) {
    stop(sprintf("expected 0 or 2 interaction features, found %d.", length(inter)), call. = FALSE)
  }

  eta
}

# internal: centre and scale to unit sd; a constant vector maps to zeros.
std_unit <- function(x) {
  s <- stats::sd(x)
  if (s == 0) rep(0, length(x)) else (x - mean(x)) / s
}
