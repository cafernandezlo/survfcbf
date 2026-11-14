#' Discretize a continuous vector into equal-frequency (quantile) bins
#'
#' Bins a continuous vector into `q` quantile (equal-frequency) bins for the FCBF
#' core, which works on discrete data. Cutpoints are learned on
#' a training fold and can be re-applied to held-out data, so the discretization
#' never uses test information.
#'
#' @details
#' With `breaks = NULL` the interior cutpoints are the `q - 1` equal-frequency
#' quantiles of `x`, learned and returned as `breaks = c(-Inf, cutpoints, Inf)`.
#' Passing those `breaks` back re-bins another vector with the same (training)
#' cutpoints; values outside the training range fall in the extreme bins.
#' Duplicate quantile cutpoints (a feature with few distinct values) are
#' collapsed, so the realised number of bins can be below `q`. The default
#' `q = max(min(floor(n / 3), 10), 2)` is Bommert et al.'s rule.
#'
#' @param x Numeric vector to discretize, without `NA`.
#' @param q Number of bins (>= 2). `NULL` (default) uses Bommert's rule from
#'   `length(x)`. Ignored when `breaks` is supplied.
#' @param breaks Optional break vector from a previous (training) call, of the
#'   form `c(-Inf, cutpoints, Inf)`, to re-apply without relearning.
#'
#' @return A list with `binned` (integer bin codes, length `length(x)`) and
#'   `breaks` (the cutpoints used, for re-application to held-out data).
#'
#' @seealso [symmetrical_uncertainty()]
#' @export
discretize_quantile <- function(x, q = NULL, breaks = NULL) {
  stopifnot("`x` must be a numeric vector without NA" = is.numeric(x) && length(x) >= 1L && !anyNA(x))

  if (is.null(breaks)) {
    n <- length(x)
    if (is.null(q)) {
      q <- max(min(floor(n / 3), 10L), 2L)
    }
    stopifnot("`q` must be a single integer >= 2" = length(q) == 1L && is.numeric(q) && is.finite(q) && q >= 2 && q == as.integer(q))
    q <- as.integer(q)
    interior <- stats::quantile(x, probs = seq_len(q - 1L) / q, names = FALSE, type = 7)
    breaks <- c(-Inf, unique(interior), Inf)        # collapse duplicate cutpoints
  } else {
    stopifnot("`breaks` must be a numeric vector of length >= 2" = is.numeric(breaks) && length(breaks) >= 2L)
  }

  binned <- as.integer(cut(x, breaks = breaks, labels = FALSE, include.lowest = TRUE))
  list(binned = binned, breaks = breaks)
}
