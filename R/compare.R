#' Bayesian signed-rank test with a region of practical equivalence
#'
#' The Bayesian counterpart to a Wilcoxon signed-rank test over paired per-dataset
#' scores (Benavoli et al. 2017), preferred over a frequentist
#' post-hoc because with ~11 datasets the latter is underpowered. Instead of a p-value it returns the posterior probability
#' that `x` is better, that `y` is better, or that the two are practically
#' equivalent, the last being the probability mass inside a region of practical
#' equivalence (ROPE) of half-width `rope`.
#'
#' @details
#' Implements the Dirichlet-process posterior over the paired differences: a
#' pseudo-observation at zero with prior weight `s` is added, weights are drawn
#' from a `Dirichlet(s, 1, ..., 1)`, and each draw is scored by the sign of the
#' pairwise means (the signed-rank statistic). Pairs where either score is `NA`
#' are dropped, so a selector that failed on a dataset does not silently count as
#' a tie. Set `higher_better = FALSE` for loss-type scores (e.g. the integrated
#' Brier score), which flips the differences so the returned `p_x` still means
#' "`x` is better".
#'
#' @section Choice of summary:
#' The posterior is a distribution *over* the three region masses, so it admits
#' two different scalar summaries, and they do not behave alike:
#'
#' * `"mean"` averages each region's mass over the draws. It is
#'   bounded above by the population quantity
#'   `P((Z_i + Z_j)/2 > rope)`, which for a real but small effect sits well below
#'   1: at a mean difference of `0.03` with `sd 0.055` the ceiling is `0.60`, so
#'   `p_x` **cannot** reach a `0.95` decision threshold no matter how many
#'   datasets are added. The ceiling is a property of this summary, not evidence
#'   that the comparison is undecidable.
#' * `"argmax"` (the default) is the convention of the Benavoli et al. reference
#'   implementation: the fraction of draws in which each region holds the largest
#'   mass. It has no such ceiling and converges to 1 for the winning region as pairs
#'   accumulate (in the same setting: `0.999` at `n = 11`, `1.000` by `n = 100`).
#'
#' Use `"argmax"` to reproduce the cited reference or to apply a `0.95` rule.
#' `"mean"` answers the different question "how much posterior mass, on average,
#' favours `x`", and it decides far fewer comparisons at the same threshold.
#'
#' @param x,y Numeric vectors of paired per-dataset scores, same length.
#' @param rope Half-width of the region of practical equivalence, in the units of
#'   the score (e.g. `0.02` on a C-index). Must be positive.
#' @param higher_better `TRUE` (default) if larger scores are better.
#' @param s Prior weight of the pseudo-observation at zero. Default `0.5`.
#' @param nsim Posterior draws. Default `50000`.
#' @param summary How the posterior over region masses is reduced to three
#'   numbers, `"argmax"` (default, the fraction of draws in which the region is
#'   largest, the Benavoli reference convention) or `"mean"` (mean mass per
#'   region). See the "Choice of summary" section, the two differ materially.
#'
#' @return A named numeric vector with `p_x`, `p_rope` and `p_y`, the summaries of
#'   the regions "`x` is better", "practically equivalent" and "`y` is better"
#'   (under `"argmax"`, the share of posterior draws in which each region holds the
#'   largest mass), and `n` (complete pairs used). All `NA` when fewer than two
#'   complete pairs.
#'
#' @references Benavoli, Corani, Demsar, Zaffalon (2017). Time for a change: a
#'   tutorial for comparing multiple classifiers through Bayesian analysis.
#'   Journal of Machine Learning Research 18(77), 1-36.
#'
#' @examples
#' set.seed(1)
#' a <- c(0.62, 0.58, 0.66, 0.61, 0.59)
#' b <- a - c(0.05, 0.03, 0.06, 0.01, 0.04)
#' bayes_signed_rank(a, b, rope = 0.02, nsim = 2000)
#'
#' @export
bayes_signed_rank <- function(x, y, rope, higher_better = TRUE, s = 0.5, nsim = 50000,
                              summary = c("argmax", "mean")) {
  summary <- match.arg(summary)
  stopifnot("`x` and `y` must be numeric vectors of the same length" =
              is.numeric(x) && is.numeric(y) && length(x) == length(y),
            "`rope` must be a single positive number" =
              length(rope) == 1L && is.numeric(rope) && is.finite(rope) && rope > 0,
            "`s` must be a single positive number" =
              length(s) == 1L && is.numeric(s) && s > 0,
            "`nsim` must be a single positive integer" =
              length(nsim) == 1L && is.numeric(nsim) && nsim >= 1)
  ok <- is.finite(x) & is.finite(y)
  z  <- (x - y)[ok]
  if (!higher_better) z <- -z                     # loss-type score: flip so p_x still means "x better"
  n <- length(z)
  if (n < 2L) return(c(p_x = NA_real_, p_rope = NA_real_, p_y = NA_real_, n = n))
  zz <- c(0, z)                                   # pseudo-observation at zero
  al <- c(s, rep(1, n))
  Mp <- outer(zz, zz, "+") / 2                    # pairwise means (the signed-rank statistic)
  Rx <- Mp >  rope
  Rr <- abs(Mp) <= rope
  Ry <- Mp < -rope
  P <- matrix(0, as.integer(nsim), 3L)
  for (b in seq_len(as.integer(nsim))) {
    g <- stats::rgamma(n + 1L, al); w <- g / sum(g)
    ww <- w %o% w
    P[b, ] <- c(sum(ww * Rx), sum(ww * Rr), sum(ww * Ry))
  }
  m <- if (summary == "mean") {
    colMeans(P)
  } else {
    tabulate(max.col(P, ties.method = "first"), nbins = 3L) / nrow(P)
  }
  c(p_x = m[1L], p_rope = m[2L], p_y = m[3L], n = n)
}
