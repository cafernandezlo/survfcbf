#' Generate proportional-hazards event times by Bender inversion
#'
#' Draws right-continuous survival times from a Weibull proportional-hazards
#' model with linear predictor `eta`, by inverting the survival function
#' (Bender, Augustin & Blettner 2005). This is the default (PH) arm of the
#' simulator.
#'
#' @details
#' With `U ~ Unif(0, 1)`, `T = (-log(U) / (lambda * exp(eta)))^(1 / a)`. Higher
#' `eta` raises the hazard and shortens `T`. The defaults `a = 1.5` and
#' `lambda = log(2)` place the baseline median (at `eta = 0`) exactly at `1`, so
#' time is measured in units of the baseline median. With `a = 1` the model
#' reduces to the exponential. Reproducible through `seed`; with `seed = NULL`
#' the current RNG state is used (so orchestration-level L'Ecuyer streams are
#' respected).
#'
#' @param eta Numeric vector of linear predictors, e.g. from
#'   [build_linear_predictor()].
#' @param shape_a Weibull shape `a` (> 0). Default `1.5`.
#' @param scale_lambda Baseline rate `lambda` (> 0). Default `log(2)`.
#' @param seed Integer or `NULL`. Seed for reproducible draws.
#'
#' @return A numeric vector of event times, the same length as `eta`.
#'
#' @seealso [aft_loglogistic_time()] for the non-proportional-hazards arm.
#' @export
bender_weibull_time <- function(eta, shape_a = 1.5, scale_lambda = log(2), seed = NULL) {
  stopifnot(
    "`eta` must be a finite numeric vector"        = is.numeric(eta) && length(eta) >= 1L && all(is.finite(eta)),
    "`shape_a` must be a single positive number"   = length(shape_a) == 1L && is.numeric(shape_a) && is.finite(shape_a) && shape_a > 0,
    "`scale_lambda` must be a single positive number" = length(scale_lambda) == 1L && is.numeric(scale_lambda) && is.finite(scale_lambda) && scale_lambda > 0
  )
  U <- with_seed(seed, stats::runif(length(eta)))
  ((-log(U)) / (scale_lambda * exp(eta)))^(1 / shape_a)
}

#' Generate non-proportional-hazards event times (log-logistic AFT)
#'
#' Draws survival times from a log-logistic accelerated-failure-time model with
#' linear predictor `eta`. This is the non-proportional-hazards arm of the
#' simulator, used to probe how the martingale-residual target
#' transform degrades when its Cox-shaped (PH) null is misspecified.
#'
#' @details
#' With `U ~ Unif(0, 1)`, the baseline `T0 = (U / (1 - U))^(1 / shape)` is
#' log-logistic with unit scale (median `1`) and `T = exp(-eta) * T0`. Higher
#' `eta` accelerates failure (shorter `T`), matching the risk direction of the
#' PH arm, but the hazard is non-monotone and the hazards of two `eta` groups are
#' not proportional, so a Cox PH model is misspecified. The baseline median (at
#' `eta = 0`) is `1`, comparable to [bender_weibull_time()].
#'
#' @param eta Numeric vector of linear predictors.
#' @param shape Log-logistic shape (> 0); values above `1` give a non-monotone
#'   hazard. Default `3`.
#' @param seed Integer or `NULL`. Seed for reproducible draws.
#'
#' @return A numeric vector of event times, the same length as `eta`.
#'
#' @seealso [bender_weibull_time()] for the proportional-hazards arm.
#' @export
aft_loglogistic_time <- function(eta, shape = 3, seed = NULL) {
  stopifnot(
    "`eta` must be a finite numeric vector"     = is.numeric(eta) && length(eta) >= 1L && all(is.finite(eta)),
    "`shape` must be a single positive number"  = length(shape) == 1L && is.numeric(shape) && is.finite(shape) && shape > 0
  )
  U <- with_seed(seed, stats::runif(length(eta)))
  T0 <- (U / (1 - U))^(1 / shape)        # log-logistic baseline, scale 1, median 1
  exp(-eta) * T0
}
