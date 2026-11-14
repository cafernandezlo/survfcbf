# internal: columns of the first distractor subblock (default disjoint covariates)
default_disjoint_cols <- function(manifest) {
  if (is.null(manifest) || !all(c("role", "cluster_id") %in% names(manifest))) {
    stop("informative_disjoint needs a manifest with role/cluster_id, or explicit `disjoint_cols`.", call. = FALSE)
  }
  d <- which(manifest$role == "distractor")
  if (length(d) == 0L) {
    stop("informative_disjoint needs distractor features in the manifest, or explicit `disjoint_cols`.", call. = FALSE)
  }
  first_id <- min(manifest$cluster_id[d], na.rm = TRUE)
  which(manifest$role == "distractor" & manifest$cluster_id == first_id)
}

# internal: per-subject censoring-rate multiplier for the exponential mechanisms
censoring_multiplier <- function(mechanism, eta, X, manifest, gamma, disjoint_cols) {
  switch(mechanism,
    independent = rep(1, length(eta)),
    informative_eta = exp(gamma * eta),
    informative_disjoint = {
      cols <- if (is.null(disjoint_cols)) default_disjoint_cols(manifest) else disjoint_cols
      exp(gamma * std_unit(rowMeans(X[, cols, drop = FALSE])))
    }
  )
}

#' Generate censoring times under independent or informative mechanisms
#'
#' Draws latent censoring times `C` for a survival simulation.
#' The observed data are then `time = pmin(T, C)` and `status = T <= C`, with `T`
#' from [bender_weibull_time()] or [aft_loglogistic_time()].
#'
#' @details
#' Two knobs are separated. The **fraction** knob is set through `param`, which
#' [calibrate_censoring()] tunes to a target censoring fraction: under
#' `fraction_model = "exponential"`, `C ~ Exp(rate = kappa * mult)`; under
#' `"administrative"`, `C = tau` (a fixed study-end cutoff, only valid for the
#' independent mechanism). The **informativeness** knob is `mechanism`:
#'
#' - `independent`: `mult = 1`, so `C` is independent of `eta` and `X` (respects
#'   the martingale / IPCW assumptions).
#' - `informative_eta`: `mult = exp(gamma * eta)`, the worst case where the
#'   censoring hazard tracks the event risk.
#' - `informative_disjoint`: `mult = exp(gamma * z)` with `z` a standardised
#'   summary of covariates disjoint from the event drivers (by default the first
#'   distractor subblock). Censoring then depends on observed covariates that do
#'   not drive the event, which biases the marginal Kaplan-Meier censoring
#'   estimator `G_hat` while conditional independence given `X` still holds.
#'
#' @param eta Numeric vector of linear predictors.
#' @param X Numeric `n` by `p` feature matrix (`nrow(X) == length(eta)`).
#' @param param Calibrated censoring parameter: `kappa` (exponential) or `tau`
#'   (administrative), typically from [calibrate_censoring()].
#' @param mechanism One of `"independent"`, `"informative_eta"`,
#'   `"informative_disjoint"`.
#' @param fraction_model `"exponential"` (default) or `"administrative"`.
#' @param manifest The [make_sigma_blocks()] manifest, used to locate the default
#'   disjoint subblock for `informative_disjoint`.
#' @param gamma Informativeness strength for the informative mechanisms. Default
#'   `0.7`.
#' @param disjoint_cols Optional integer columns of `X` to drive
#'   `informative_disjoint`, overriding the default distractor subblock.
#' @param seed Integer or `NULL`. Seed for reproducible draws.
#'
#' @return A numeric vector of censoring times, the same length as `eta`.
#'
#' @seealso [calibrate_censoring()]
#' @export
make_censoring <- function(eta, X, param,
                           mechanism = c("independent", "informative_eta", "informative_disjoint"),
                           fraction_model = c("exponential", "administrative"),
                           manifest = NULL, gamma = 0.7, disjoint_cols = NULL, seed = NULL) {
  mechanism <- match.arg(mechanism)
  fraction_model <- match.arg(fraction_model)
  stopifnot(
    "`eta` must be a finite numeric vector"                  = is.numeric(eta) && length(eta) >= 1L && all(is.finite(eta)),
    "`X` must be a numeric matrix with nrow == length(eta)"  = is.matrix(X) && is.numeric(X) && nrow(X) == length(eta),
    "`param` must be a single positive number"               = length(param) == 1L && is.numeric(param) && is.finite(param) && param > 0,
    "`gamma` must be a single finite number"                 = length(gamma) == 1L && is.numeric(gamma) && is.finite(gamma)
  )
  n <- length(eta)

  if (fraction_model == "administrative") {
    if (mechanism != "independent") {
      stop("administrative censoring supports only mechanism = 'independent'.", call. = FALSE)
    }
    return(rep(param, n))                       # C = tau (administrative cutoff)
  }

  mult <- censoring_multiplier(mechanism, eta, X, manifest, gamma, disjoint_cols)
  E_c <- with_seed(seed, stats::rexp(n))
  E_c / (param * mult)                          # C ~ Exp(rate = kappa * mult)
}

#' Calibrate the censoring parameter to a target censoring fraction
#'
#' Finds the censoring parameter (`kappa` for the exponential model, `tau` for
#' the administrative one) at which the realised censoring fraction
#' `mean(C < T)` hits `target_fraction`, given an event-time
#' pre-sample.
#'
#' @details
#' For the exponential model the censoring base draws are fixed once, so
#' `fraction(kappa) = mean(E_c < kappa * mult * T)` is monotone increasing and
#' deterministic in `kappa`; the bracket is seeded with the closed-form anchor
#' `kappa = log(2) * f / (1 - f)` (exact when `T ~ Exp(log 2)`) and refined by
#' bisection. The administrative `tau` is the `(1 - f)` quantile of the event
#' times directly, with no bisection.
#'
#' @param eta,X,manifest,mechanism,fraction_model,gamma,disjoint_cols As in
#'   [make_censoring()]; the mechanism and covariates must match those used to
#'   later generate the censoring.
#' @param event_times Numeric vector of event times `T` (same length as `eta`),
#'   e.g. from [bender_weibull_time()].
#' @param target_fraction Target censoring fraction, in `(0, 1)`.
#' @param tol Tolerance on the fraction for the bisection. Default `0.005`.
#' @param max_iter Maximum bisection iterations. Default `80`.
#' @param seed Integer or `NULL`. Seed for the reproducible censoring base draws.
#'
#' @return The calibrated parameter (`kappa` or `tau`), a single numeric.
#'
#' @seealso [make_censoring()]
#' @export
calibrate_censoring <- function(eta, X, event_times, target_fraction,
                                mechanism = c("independent", "informative_eta", "informative_disjoint"),
                                fraction_model = c("exponential", "administrative"),
                                manifest = NULL, gamma = 0.7, disjoint_cols = NULL,
                                tol = 0.005, max_iter = 80, seed = NULL) {
  mechanism <- match.arg(mechanism)
  fraction_model <- match.arg(fraction_model)
  stopifnot(
    "`target_fraction` must be a single number in (0, 1)" = length(target_fraction) == 1L && is.numeric(target_fraction) && target_fraction > 0 && target_fraction < 1,
    "`event_times` must be numeric with length(eta)"      = is.numeric(event_times) && length(event_times) == length(eta),
    "`tol` must be a single positive number"              = length(tol) == 1L && is.numeric(tol) && tol > 0
  )

  if (fraction_model == "administrative") {
    if (mechanism != "independent") {
      stop("administrative censoring supports only mechanism = 'independent'.", call. = FALSE)
    }
    return(as.numeric(stats::quantile(event_times, probs = 1 - target_fraction, names = FALSE)))
  }

  mult <- censoring_multiplier(mechanism, eta, X, manifest, gamma, disjoint_cols)
  E_c <- with_seed(seed, stats::rexp(length(eta)))
  fraction_at <- function(kappa) mean(E_c < kappa * mult * event_times)

  ## bracket: closed-form anchor for the clean exponential case, then expand
  hi <- log(2) * target_fraction / (1 - target_fraction)
  if (!is.finite(hi) || hi <= 0) hi <- 1
  expand <- 0L
  while (fraction_at(hi) < target_fraction && expand < 80L) {
    hi <- hi * 2
    expand <- expand + 1L
  }

  ## bisection on [0, hi]
  lo <- 0
  mid <- hi
  for (it in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    f <- fraction_at(mid)
    if (abs(f - target_fraction) < tol) break
    if (f < target_fraction) lo <- mid else hi <- mid
  }
  mid
}
