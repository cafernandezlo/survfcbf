#' Transform a right-censored outcome into an uncensored pseudo-target
#'
#' Maps `(time, status)` to a continuous scalar that the FCBF relevance leg can
#' use, the standard route for feeding mutual-information filters survival data.
#' The default is the martingale residual of a covariate-free
#' Cox model, which is the transform Bommert et al. (2022) use and that this
#' package reuses rather than claims as a contribution.
#'
#' @details
#' For `method = "martingale"` a null Cox model `coxph(Surv(time, status) ~ 1)`
#' is fitted and its martingale residuals returned: `status - H0(time)`, where
#' `H0` is the baseline (Nelson-Aalen) cumulative hazard. They sum to zero, an
#' early event gives a residual near `+1` while a long-surviving subject gives a
#' negative one, and two subjects observed at the same time differ by exactly
#' their status. Symmetrical uncertainty is direction-agnostic, so only the
#' ordering matters downstream.
#'
#' When fitted inside a cross-validation fold the null Cox model must be fit on
#' the training fold only; [survfcbf_select()] enforces that when it runs inside
#' [nested_cv_select()].
#'
#' @param time Numeric vector of follow-up times (all `> 0`).
#' @param status Event indicator (`1` = event, `0` = censored; logical accepted),
#'   the same length as `time`.
#' @param method Target transform. `"martingale"` (default) returns the null-Cox
#'   martingale residuals; `"ipcw_binned"` returns the observed `time` as the
#'   target, to be discretized and IPCW-weighted by [survfcbf_select()] (the 5.5
#'   ablation).
#'
#' @return A numeric pseudo-target, the same length as `time`.
#'
#' @export
survival_target <- function(time, status, method = c("martingale", "ipcw_binned")) {
  method <- match.arg(method)
  if (is.logical(status)) status <- as.integer(status)
  stopifnot(
    "`time` must be a numeric vector (length >= 2) of positive values without NA" =
      is.numeric(time) && length(time) >= 2L && !anyNA(time) && all(time > 0),
    "`status` must be 0/1, the same length as `time`, without NA" =
      length(status) == length(time) && !anyNA(status) && all(status %in% c(0L, 1L))
  )

  if (method == "ipcw_binned") {
    return(as.numeric(time))     # observed time as the target; IPCW weighting is applied downstream
  }

  fit <- survival::coxph(survival::Surv(time, status) ~ 1)
  as.numeric(stats::residuals(fit, type = "martingale"))
}
