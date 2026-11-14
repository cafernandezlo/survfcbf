#' Sample correlated features with chosen marginals via a Gaussian copula
#'
#' Draws `n` observations of the `p` latent-correlated features encoded by
#' `Sigma` and maps them to the requested marginal distribution with a Gaussian
#' copula, returning the realised feature matrix and (optionally) its empirical
#' correlation.
#'
#' @details
#' Latent draws `Z ~ N(0, Sigma)` exploit the block-diagonal structure of
#' `Sigma`: each correlated block is sampled through its own Cholesky factor and
#' the independent columns (the identity-covariance noise and interaction) by
#' plain standard normals, avoiding a dense `O(p^3)` decomposition. For
#' `marginals = "gaussian"` (a sanity mode) the features are `X = Z`. Otherwise
#' the Gaussian copula `U = pnorm(Z)` is pushed through a per-feature inverse
#' cdf: a negative-binomial quantile (`marginals = "negbin"`) or the empirical
#' quantile of a real reference gene (`marginals = "empirical"`, using a
#' reference from [build_empirical_reference()]).
#'
#' Because the copula maps each margin through a monotone transform, the rank
#' (and hence quantile-binned) structure of the latent block correlations is
#' preserved, while the realised linear correlation of non-Gaussian margins is
#' attenuated relative to the latent `rho`. The realised correlation is returned
#' so this attenuation can be measured.
#'
#' Sampling reproducibility: with `seed = NULL` (default) features are drawn from
#' the current RNG state, so per-replicate L'Ecuyer-CMRG streams set up by the
#' orchestration layer are respected; pass an integer `seed` for standalone
#' reproducibility (the global RNG state is restored afterwards).
#'
#' @param n Integer. Number of observations (rows of `X`).
#' @param Sigma Numeric `p` by `p` positive-definite latent covariance, e.g. the
#'   `Sigma` from [make_sigma_blocks()].
#' @param marginals Marginal family: `"gaussian"`, `"negbin"` or `"empirical"`.
#' @param negbin_mu,negbin_size Mean and size (dispersion) of the negative
#'   binomial margin, used when `marginals = "negbin"`. Defaults `50` and `2`.
#' @param reference A `survfcbf_empirical_reference` from
#'   [build_empirical_reference()] whose `p` matches `ncol(Sigma)`, used when
#'   `marginals = "empirical"`.
#' @param seed Integer or `NULL`. Seed for reproducible sampling.
#' @param return_cor Logical. If `TRUE` (default) include the realised `p` by `p`
#'   Pearson correlation `cor(X)`; set `FALSE` to skip it for large `p`.
#'
#' @return A list with `X` (an `n` by `p` numeric matrix) and `realized_cor`
#'   (the `p` by `p` Pearson correlation of `X`, or `NULL` when `return_cor` is
#'   `FALSE`).
#'
#' @seealso [make_sigma_blocks()], [build_empirical_reference()]
#' @export
sample_features <- function(n, Sigma,
                            marginals = c("gaussian", "negbin", "empirical"),
                            negbin_mu = 50, negbin_size = 2,
                            reference = NULL, seed = NULL, return_cor = TRUE) {
  marginals <- match.arg(marginals)
  stopifnot(
    "`n` must be a single positive integer"    = length(n) == 1L && is.numeric(n) && is.finite(n) && n == as.integer(n) && n >= 1,
    "`Sigma` must be a square numeric matrix"  = is.matrix(Sigma) && is.numeric(Sigma) && nrow(Sigma) == ncol(Sigma)
  )
  n <- as.integer(n)
  p <- ncol(Sigma)

  if (marginals == "empirical") {
    if (is.null(reference) || !inherits(reference, "survfcbf_empirical_reference")) {
      stop("`marginals = \"empirical\"` needs a `reference` from build_empirical_reference().",
           call. = FALSE)
    }
    if (reference$p != p) {
      stop(sprintf("reference is built for p = %d but Sigma has p = %d.", reference$p, p),
           call. = FALSE)
    }
  }

  blocks <- sigma_blocks(Sigma)
  Z <- with_seed(seed, sample_block_normal(n, Sigma, blocks))

  X <- switch(marginals,
    gaussian = Z,
    negbin = {
      U <- stats::pnorm(Z)
      matrix(stats::qnbinom(U, mu = negbin_mu, size = negbin_size), nrow = n, ncol = p)
    },
    empirical = {
      U <- stats::pnorm(Z)
      empirical_quantile_matrix(U, reference$quantiles)
    }
  )
  colnames(X) <- colnames(Sigma)

  realized_cor <- if (isTRUE(return_cor)) stats::cor(X) else NULL
  list(X = X, realized_cor = realized_cor)
}

# internal: connected components (blocks) of Sigma's non-zero correlation graph,
# via union-find over the off-diagonal non-zeros
sigma_blocks <- function(Sigma) {
  p <- ncol(Sigma)
  parent <- seq_len(p)
  root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  nz <- which(Sigma != 0, arr.ind = TRUE)
  nz <- nz[nz[, 1L] < nz[, 2L], , drop = FALSE]        # off-diagonal upper pairs
  for (k in seq_len(nrow(nz))) {
    a <- root(nz[k, 1L])
    b <- root(nz[k, 2L])
    if (a != b) parent[a] <- b
  }
  ids <- vapply(seq_len(p), root, integer(1))
  unname(split(seq_len(p), ids))
}

# internal: sample n rows of N(0, Sigma) block by block -- a Cholesky factor per
# correlated block, scaled standard normals for singletons -- exploiting the
# block-diagonal structure to avoid a dense p x p decomposition
sample_block_normal <- function(n, Sigma, blocks) {
  p <- ncol(Sigma)
  Z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  for (b in blocks) {
    if (length(b) == 1L) {
      v <- Sigma[b, b]
      if (v != 1) Z[, b] <- Z[, b] * sqrt(v)
    } else {
      Z[, b] <- Z[, b, drop = FALSE] %*% chol(Sigma[b, b])
    }
  }
  Z
}
