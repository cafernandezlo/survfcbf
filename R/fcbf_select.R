#' FCBF predominant-feature selection on discrete data
#'
#' The two-stage Fast Correlation-Based Filter (Yu and Liu, 2003) on a discrete
#' feature matrix and discrete target. A relevance stage keeps
#' features whose symmetrical uncertainty with the target reaches `delta` and
#' ranks them; a redundancy stage prunes by an approximate Markov blanket. The
#' survivors are the predominant features. The result carries a full trace of the
#' Markov-blanket decisions, which is the instrumentation the method is built to
#' expose.
#'
#' @details
#' Relevance: keep features with `SU(Fi, y) >= delta`, sorted descending by that
#' relevance. Redundancy: stepping down the ranked list, a more relevant `Fi`
#' removes a less relevant `Fj` when `SU(Fi, Fj) >= SU(Fj, y)`, i.e. `Fi` is an
#' approximate Markov blanket for `Fj`. Pairwise `SU(Fi, Fj)` is computed lazily,
#' only for pairs the pruning actually reaches.
#'
#' A positive `top` switches the cardinality control: instead of letting `delta`
#' decide how many features survive, it walks down the relevance ranking keeping
#' each feature that no already-kept (more relevant) feature blankets, and stops
#' once `top` survive. This returns the `top` most relevant predominant features
#' FCBF would select (the early stop is exact, since a feature's survival depends
#' only on more relevant ones), with the same SU Markov-blanket pruning. It is the
#' mode for the low-event / high-dimensional regime, where the permutation-null
#' `delta` (the max-SU upper tail over many features) inflates above the signal and
#' the relevance stage would otherwise select nothing.
#'
#' The `>=` follows Yu and Liu (2003); the Bioconductor `FCBF` package uses a
#' strict `>`, so the two selections diverge only at exact SU ties (for instance
#' a feature identical to the target), where `>=` prunes the duplicate and `>`
#' keeps it. On non-degenerate data the selections are identical.
#'
#' Inputs must already be discrete (the survival target transform and the
#' quantile discretisation happen upstream, see [survfcbf_select()]). For efficiency on wide data use a
#' positive `delta`, so the redundancy pass runs over the relevant set rather
#' than all columns.
#'
#' @param X A discrete feature matrix or data.frame, observations in rows,
#'   features in columns.
#' @param y A discrete target vector, length `nrow(X)`. Required unless
#'   `relevance` is supplied (then it is unused, as relevance is precomputed).
#' @param delta Relevance threshold, in `[0, 1]`. Default `0` (no relevance
#'   filter, pure redundancy pruning).
#' @param relevance Optional precomputed per-feature relevance (length
#'   `ncol(X)`, on the same `[0, 1]` scale as SU), e.g. a rescaled univariate
#'   Cox concordance. `NULL` (default) uses `SU(Fi, y)`. The redundancy stage
#'   always uses SU between features, so the Markov-blanket inequality stays
#'   SU-comparable.
#' @param weights Optional non-negative observation weights for the SU
#'   computations (e.g. IPCW weights), length `nrow(X)`. `NULL` (default) is
#'   unweighted.
#' @param top Optional positive integer cardinality cap. When supplied, the
#'   redundancy stage stops once `top` predominant features survive and `delta`
#'   acts only as an optional relevance floor (default `0`, no floor). `NULL`
#'   (default) prunes the whole relevant set, classic FCBF with the set size
#'   determined by `delta`.
#' @param tiebreak Optional numeric secondary key, length `ncol(X)`, that orders
#'   features with equal relevance (higher first). With coarse discretisation many
#'   features share the same SU, and without a secondary key their order, and so
#'   the selection, depends on the column order of `X`. `NULL` (default) keeps the
#'   classic behaviour of breaking ties by column position, as the Bioconductor
#'   `FCBF` package does. [survfcbf_select()] supplies the absolute Spearman
#'   correlation of each continuous feature with the survival pseudo-target.
#'
#' @return A list with `selected` (integer column indices of the predominant
#'   features, relevance-ordered) and `trace`, itself a list of:
#'   \describe{
#'     \item{`relevance_su`}{named numeric `SU(Fi, y)` for every feature.}
#'     \item{`removed`}{a data.frame with one row per pruned feature: `feature`,
#'       the `blanket` that removed it, and the two SU values
#'       `su_blanket_feature` (`SU(blanket, feature)`) and `su_feature_target`
#'       (`SU(feature, y)`).}
#'     \item{`representative_per_cluster`}{a named list mapping each predominant
#'       feature to the group it represents (itself plus the features it
#'       blanketed).}
#'   }
#'   Under a positive `top` the `removed` and `representative_per_cluster` traces
#'   cover only the features examined before the cardinality cap was reached;
#'   features further down the ranking are neither selected nor traced.
#'
#' @seealso [symmetrical_uncertainty()]
#'
#' @examples
#' set.seed(1)
#' y  <- sample(1:3, 300, replace = TRUE)
#' f1 <- ifelse(runif(300) < 0.1, sample(1:3, 300, replace = TRUE), y)   # relevant
#' f2 <- ifelse(runif(300) < 0.1, sample(1:3, 300, replace = TRUE), f1)  # near-copy of f1
#' f3 <- sample(1:3, 300, replace = TRUE)                                 # irrelevant
#' X  <- cbind(f1 = f1, f2 = f2, f3 = f3)
#'
#' res <- fcbf_select(X, y, delta = 0.05)
#' res$selected        # f1 only
#' res$trace$removed   # f2 pruned, with f1 as its approximate Markov blanket
#' @export
fcbf_select <- function(X, y = NULL, delta = 0, relevance = NULL, weights = NULL, top = NULL,
                        tiebreak = NULL) {
  stopifnot(
    "`X` must be a matrix or data.frame"        = is.matrix(X) || is.data.frame(X),
    "`delta` must be a single number in [0, 1]" = length(delta) == 1L && is.numeric(delta) && delta >= 0 && delta <= 1,
    "`top` must be NULL or a positive integer"  = is.null(top) || (length(top) == 1L && is.numeric(top) && !is.na(top) && top >= 1),
    "`tiebreak` must be NULL or a numeric vector of length ncol(X) without NA" =
      is.null(tiebreak) || (is.numeric(tiebreak) && length(tiebreak) == ncol(X) && !anyNA(tiebreak))
  )
  if (!is.null(top)) top <- as.integer(top)
  p <- ncol(X)
  get_col <- function(j) X[, j]
  feat_names <- if (!is.null(colnames(X))) colnames(X) else as.character(seq_len(p))

  ## relevance stage (precomputed `relevance` overrides the default SU(Fi, y))
  if (is.null(relevance)) {
    stopifnot("`y` must have length nrow(X) when `relevance` is not supplied" = !is.null(y) && length(y) == nrow(X))
    su_y <- vapply(seq_len(p), function(j) symmetrical_uncertainty(get_col(j), y, weights), numeric(1))
  } else {
    stopifnot("`relevance` must be numeric of length ncol(X), in [0, 1]" =
                is.numeric(relevance) && length(relevance) == p &&
                !anyNA(relevance) && all(relevance >= 0 & relevance <= 1))
    su_y <- relevance
  }
  ord <- which(su_y >= delta)
  ## relevance-ordered feature indices. order() is stable, so features tied on both keys
  ## keep their column order.
  key2 <- if (is.null(tiebreak)) numeric(p) else tiebreak
  ord <- ord[order(su_y[ord], key2[ord], decreasing = TRUE)]

  ## redundancy stage: approximate Markov-blanket pruning down the ranked list
  n_ord <- length(ord)
  rm_feature <- integer(0)
  rm_blanket <- integer(0)
  rm_su_pq <- numeric(0)
  rm_su_qy <- numeric(0)

  if (is.null(top)) {
    ## classic FCBF: prune the whole relevant set, cardinality set by `delta`
    alive <- rep(TRUE, n_ord)
    for (a in seq_len(n_ord)) {
      if (!alive[a]) next
      Fp <- ord[a]
      b <- a + 1L
      while (b <= n_ord) {
        if (alive[b]) {
          Fq <- ord[b]
          su_pq <- symmetrical_uncertainty(get_col(Fp), get_col(Fq), weights)
          if (su_pq >= su_y[Fq]) {                      # Fp is an approx Markov blanket for Fq
            alive[b] <- FALSE
            rm_feature <- c(rm_feature, Fq)
            rm_blanket <- c(rm_blanket, Fp)
            rm_su_pq <- c(rm_su_pq, su_pq)
            rm_su_qy <- c(rm_su_qy, su_y[Fq])
          }
        }
        b <- b + 1L
      }
    }
    selected <- ord[alive]
  } else {
    ## explicit-cardinality walk: keep predominant features until `top` survive.
    ## A candidate Fq is kept unless an already-kept (more relevant) feature is an
    ## approximate Markov blanket for it; the first such blanket (the most relevant,
    ## as `kept` is relevance-ordered) is recorded, matching the classic trace.
    kept <- integer(0)                                  # positions in `ord` of survivors
    for (a in seq_len(n_ord)) {
      if (length(kept) >= top) break
      Fq <- ord[a]
      blanket <- NA_integer_
      su_blanket <- NA_real_
      for (kp in kept) {
        Fp <- ord[kp]
        su_pq <- symmetrical_uncertainty(get_col(Fp), get_col(Fq), weights)
        if (su_pq >= su_y[Fq]) { blanket <- Fp; su_blanket <- su_pq; break }
      }
      if (is.na(blanket)) {
        kept <- c(kept, a)
      } else {
        rm_feature <- c(rm_feature, Fq)
        rm_blanket <- c(rm_blanket, blanket)
        rm_su_pq <- c(rm_su_pq, su_blanket)
        rm_su_qy <- c(rm_su_qy, su_y[Fq])
      }
    }
    selected <- ord[kept]
  }

  representative_per_cluster <- stats::setNames(
    lapply(selected, function(fp) c(fp, rm_feature[rm_blanket == fp])),
    as.character(selected)
  )

  trace <- list(
    relevance_su = stats::setNames(su_y, feat_names),
    removed = data.frame(
      feature = rm_feature,
      blanket = rm_blanket,
      su_blanket_feature = rm_su_pq,
      su_feature_target = rm_su_qy
    ),
    representative_per_cluster = representative_per_cluster
  )

  list(selected = selected, trace = trace)
}
