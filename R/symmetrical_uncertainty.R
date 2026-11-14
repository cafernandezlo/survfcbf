#' Shannon entropy of a discrete vector
#'
#' @param x A discrete vector (factor, integer or character) without `NA`; its
#'   distinct values are the categories.
#' @param w Optional non-negative observation weights (length `length(x)`),
#'   for example IPCW weights. `NULL` (default) weights every observation
#'   equally.
#' @return The entropy in bits.
#' @seealso [conditional_entropy()], [symmetrical_uncertainty()]
#' @export
entropy <- function(x, w = NULL) {
  stopifnot("`x` must be a non-empty vector without NA" = length(x) >= 1L && !anyNA(x))
  if (is.null(w)) {
    px <- table(x) / length(x)
  } else {
    stopifnot("`w` must be non-negative weights of length(x) with positive sum" =
                length(w) == length(x) && !anyNA(w) && all(w >= 0) && sum(w) > 0)
    px <- tapply(w, x, sum)
    px <- px / sum(px)
  }
  px <- px[px > 0]
  -sum(px * log2(px))
}

#' Conditional entropy of two discrete vectors
#'
#' Computes `H(X | Y) = H(X, Y) - H(Y)` for discrete `x` and `y`.
#'
#' @param x,y Discrete vectors of equal length without `NA`.
#' @param w Optional non-negative observation weights (length `length(x)`).
#'   `NULL` (default) weights every observation equally.
#' @return `H(X | Y)` in bits.
#' @seealso [entropy()], [symmetrical_uncertainty()]
#' @export
conditional_entropy <- function(x, y, w = NULL) {
  stopifnot(
    "`x` and `y` must be equal-length vectors without NA" =
      length(x) == length(y) && length(x) >= 1L && !anyNA(x) && !anyNA(y)
  )
  if (is.null(w)) {
    n <- length(x)
    joint <- table(x, y) / n
    py <- table(y) / n
  } else {
    stopifnot("`w` must be non-negative weights of length(x) with positive sum" =
                length(w) == length(x) && !anyNA(w) && all(w >= 0) && sum(w) > 0)
    sw <- sum(w)
    joint <- tapply(w, list(x, y), sum)
    joint[is.na(joint)] <- 0
    joint <- joint / sw
    py <- tapply(w, y, sum) / sw
  }
  joint <- joint[joint > 0]
  py <- py[py > 0]
  (-sum(joint * log2(joint))) - (-sum(py * log2(py)))     # H(X, Y) - H(Y)
}

#' Symmetrical uncertainty of two discrete vectors
#'
#' The entropy-based, range-normalised association used by FCBF:
#' `SU(X, Y) = 2 * IG(X | Y) / (H(X) + H(Y))`, where the information gain
#' `IG(X | Y) = H(X) - H(X | Y)` is the mutual information. SU is symmetric and
#' lies in `[0, 1]` (1 for a deterministic relationship, 0 for independence); the
#' logarithm base cancels in the ratio. Two constant vectors return 0.
#'
#' @param x,y Discrete vectors of equal length without `NA`.
#' @param w Optional non-negative observation weights (length `length(x)`), e.g.
#'   IPCW weights for the censored-data ablation. `NULL` (default) is unweighted.
#' @return The symmetrical uncertainty, a single number in `[0, 1]`.
#' @seealso [entropy()], [conditional_entropy()]
#' @export
symmetrical_uncertainty <- function(x, y, w = NULL) {
  hx <- entropy(x, w)
  hy <- entropy(y, w)
  denom <- hx + hy
  if (denom == 0) {
    return(0)                                  # both constant: no shared information
  }
  ig <- hx - conditional_entropy(x, y, w)      # mutual information
  min(max(2 * ig / denom, 0), 1)               # clamp tiny numerical excursions to [0, 1]
}
