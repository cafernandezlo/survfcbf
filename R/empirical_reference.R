#' Build a per-feature empirical marginal reference for the simulator
#'
#' Turns a real expression matrix (genes by samples) into a reusable reference
#' of `p` per-feature inverse cumulative distribution functions, so
#' [sample_features()] can impose realistic, gene-like marginals on the
#' latent-correlated features via a Gaussian copula. The reference preparation
#' (gene filtering and gene-to-feature assignment) happens once here and is
#' reused across simulator replicates.
#'
#' @details
#' Genes enter a usable pool if they have positive variance and are expressed
#' above `min_expr` in at least `min_frac` of samples. `p` pool genes are then
#' drawn (with a fixed `seed` for reproducibility) and assigned one-to-one to the
#' simulated features; each assigned gene contributes its observed sample values
#' as that feature's empirical reference. The assignment is independent of any
#' feature role, so marginal variance is not confounded with relevance. If `p`
#' exceeds the pool size the genes are drawn with replacement (with a warning).
#'
#' The `reference` matrix must already be on the desired marginal scale (for
#' example `log2(TPM + 1)`); no transformation is applied here. Extracting that
#' matrix from a `SummarizedExperiment` is intentionally kept out of the package
#' (see `data-raw/prepare_tcga_reference.R`) so the package depends only on base
#' R.
#'
#' @param reference Numeric matrix, genes (rows) by samples (columns), on the
#'   target marginal scale. Row names, if present, are kept as reference gene ids.
#' @param p Integer. Number of features to build references for (the `p` of the
#'   `Sigma` later passed to [sample_features()]).
#' @param seed Integer or `NULL`. Seed for the reproducible gene-to-feature
#'   assignment.
#' @param min_expr Numeric. A sample counts a gene as expressed when its value
#'   exceeds this threshold (on the reference scale). Default `0`.
#' @param min_frac Numeric in `[0, 1]`. Minimum fraction of samples in which a
#'   gene must be expressed for it to enter the usable pool. Default `0.5`.
#' @param low_var_features Optional integer feature indices (typically the driver
#'   cluster columns) to receive the lowest-variance usable genes, for the
#'   low-variance-signal scenario; the remaining features keep the random
#'   assignment. `NULL` (default) assigns every feature at random.
#'
#' @return An object of class `survfcbf_empirical_reference`: a list with
#'   `quantiles` (a samples by `p` numeric matrix whose column `j` holds the
#'   sorted reference values for feature `j`), `genes` (the assigned reference
#'   gene ids, length `p`), `p`, `n_samples`, and `n_pool` (usable pool size).
#'
#' @seealso [sample_features()]
#' @export
build_empirical_reference <- function(reference, p, seed = NULL,
                                      min_expr = 0, min_frac = 0.5,
                                      low_var_features = NULL) {
  stopifnot(
    "`reference` must be a numeric matrix"          = is.matrix(reference) && is.numeric(reference),
    "`p` must be a single positive integer"         = length(p) == 1L && is.numeric(p) && is.finite(p) && p == as.integer(p) && p >= 1,
    "`min_expr` must be a single finite number"     = length(min_expr) == 1L && is.numeric(min_expr) && is.finite(min_expr),
    "`min_frac` must be a single number in [0, 1]"  = length(min_frac) == 1L && is.numeric(min_frac) && min_frac >= 0 && min_frac <= 1
  )
  p <- as.integer(p)
  n_samples <- ncol(reference)
  if (n_samples < 2L) {
    stop("`reference` needs at least 2 samples (columns).", call. = FALSE)
  }

  ## usable gene pool: positive variance and expressed in >= min_frac of samples
  expressed_frac <- rowMeans(reference > min_expr, na.rm = TRUE)
  row_var <- apply(reference, 1L, stats::var, na.rm = TRUE)
  usable <- which(is.finite(row_var) & row_var > 0 & expressed_frac >= min_frac)
  n_pool <- length(usable)
  if (n_pool == 0L) {
    stop("No usable genes after filtering; relax `min_expr` / `min_frac`.", call. = FALSE)
  }

  ## assign p genes to features (fixed seed), with replacement only if needed
  use_replace <- p > n_pool
  if (use_replace) {
    warning(sprintf("p = %d exceeds the usable gene pool (%d); sampling genes with replacement.",
                    p, n_pool), call. = FALSE)
  }
  pick <- with_seed(seed, {
    if (is.null(low_var_features)) {
      sample(usable, size = p, replace = use_replace)
    } else {
      stopifnot(
        "`low_var_features` must be feature indices in 1..p" =
          is.numeric(low_var_features) && all(low_var_features >= 1 & low_var_features <= p)
      )
      lv <- unique(as.integer(low_var_features))
      low_genes <- usable[order(row_var[usable])][seq_len(min(length(lv), n_pool))]
      lv <- lv[seq_along(low_genes)]
      rest_pool <- setdiff(usable, low_genes)            # other features avoid the low-var genes
      if (length(rest_pool) == 0L) {                     # low-var assignment consumed the whole pool
        warning("the low-variance assignment used every usable gene; remaining features ",
                "reuse genes with replacement.", call. = FALSE)
        rest_pool <- usable
      }
      rest_feats <- setdiff(seq_len(p), lv)
      chosen <- integer(p)
      chosen[lv] <- low_genes                            # lowest-variance genes to the named features
      if (length(rest_feats) > 0L) {
        chosen[rest_feats] <- sample(rest_pool, length(rest_feats),
                                     replace = length(rest_feats) > length(rest_pool))
      }
      chosen
    }
  })

  ## sorted per-feature reference columns (samples x p) for fast quantile inversion.
  ## apply(., 1, sort) returns a samples x p matrix (one sorted column per feature),
  ## or a length-samples vector when p == 1.
  sorted_cols <- apply(reference[pick, , drop = FALSE], 1L, sort)
  qmat <- if (is.null(dim(sorted_cols))) matrix(sorted_cols, ncol = 1L) else sorted_cols

  gene_ids <- rownames(reference)
  gene_ids <- if (is.null(gene_ids)) paste0("gene", pick) else gene_ids[pick]

  structure(
    list(quantiles = qmat, genes = gene_ids, p = p,
         n_samples = n_samples, n_pool = n_pool),
    class = "survfcbf_empirical_reference"
  )
}

#' @export
print.survfcbf_empirical_reference <- function(x, ...) {
  cat(sprintf("<survfcbf empirical reference> %d features from a pool of %d genes, %d samples each\n",
              x$p, x$n_pool, x$n_samples))
  invisible(x)
}
