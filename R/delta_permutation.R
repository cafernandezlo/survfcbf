#' Relevance threshold by a permutation null of the target
#'
#' Chooses the FCBF relevance threshold `delta` from a permutation null:
#' the `prob` quantile of the maximum `SU(feature, permuted target)`
#' over `B` permutations of the discretized target. Selecting features with
#' `SU(Fi, y) >= delta` then controls the family-wise false-selection rate near
#' `1 - prob` under the null.
#'
#' @details
#' A permutation preserves each feature's and the target's marginal, so
#' `H(X_j)` and `H(y)` are constant across permutations and are precomputed once;
#' each permutation only re-evaluates the joint entropies. Fit on the training
#' fold only and seeded for reproducibility. The alternative to the permutation
#' null is a fixed `delta` grid tuned by inner cross-validation, which is simply
#' passed to [fcbf_select()] / `survfcbf_select` directly rather than estimated
#' here.
#'
#' @param X Discrete feature matrix (n by p), e.g. from [discretize_quantile()].
#' @param y Discrete target of length `nrow(X)` (the discretized pseudo-target).
#' @param B Number of target permutations. Default `200`.
#' @param prob Quantile of the per-permutation maximum SU to use as `delta`.
#'   Default `0.95`.
#' @param seed Integer or `NULL`. Seed for the reproducible permutations.
#' @param weights Optional non-negative observation weights for the SU
#'   computations (e.g. IPCW weights), length `nrow(X)`. `NULL` (default) is
#'   unweighted.
#'
#' @return The relevance threshold `delta`, a single number in `[0, 1]`.
#'
#' @seealso [fcbf_select()], [symmetrical_uncertainty()]
#' @export
delta_permutation <- function(X, y, B = 200L, prob = 0.95, seed = NULL, weights = NULL) {
  stopifnot(
    "`X` must be a matrix"                          = is.matrix(X) && nrow(X) >= 2L,
    "`y` must have length nrow(X)"                  = length(y) == nrow(X),
    "`B` must be a single positive integer"         = length(B) == 1L && is.numeric(B) && is.finite(B) && B >= 1 && B == as.integer(B),
    "`prob` must be a single number in (0, 1)"      = length(prob) == 1L && is.numeric(prob) && prob > 0 && prob < 1
  )
  B <- as.integer(B)
  p <- ncol(X)

  Hx <- apply(X, 2L, entropy, w = weights)   # invariant under target permutation

  max_su <- with_seed(seed, vapply(seq_len(B), function(b) {
    yp <- sample(y)
    Hy <- entropy(yp, weights)               # recomputed per permutation: the weighted H(yp)
    denom <- Hx + Hy                          # varies when weights are tied to subjects, not labels
    su <- vapply(seq_len(p), function(j) {
      if (denom[j] <= 0) return(0)
      min(max(2 * (Hx[j] - conditional_entropy(X[, j], yp, weights)) / denom[j], 0), 1)
    }, numeric(1))
    max(su)
  }, numeric(1)))

  stats::quantile(max_su, probs = prob, names = FALSE)
}
