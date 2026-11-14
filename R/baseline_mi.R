# internal: survfcbf's exact relevance substrate for the direct MI family. Fits the
# martingale pseudo-target and the equal-frequency cutpoints on this (training) fold,
# then returns the binned features and target as factors so praznik treats them as
# categorical and does NOT re-discretize. Sharing this with survfcbf_select is the
# whole point: only the selection rule (MI ranking vs FCBF predominance) then differs.
mi_substrate <- function(time, status, X, q = NULL) {
  stopifnot("`X` must be a numeric matrix" = is.matrix(X) && is.numeric(X),
            "`X`, `time`, `status` must share nrow/length" =
              nrow(X) == length(time) && length(time) == length(status))
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(q)) q <- max(min(floor(n / 3), 10L), 2L)
  ym <- survival_target(time, status, method = "martingale")
  y_disc <- discretize_quantile(ym, q = q)$binned
  X_disc <- vapply(seq_len(p), function(j) discretize_quantile(X[, j], q = q)$binned, integer(n))
  Xf <- as.data.frame(lapply(seq_len(p), function(j) factor(X_disc[, j])))
  names(Xf) <- if (!is.null(colnames(X))) colnames(X) else as.character(seq_len(p))
  list(Xf = Xf, yf = factor(y_disc), p = p, q = q)
}

# internal: shared praznik wrapper for the MI family. `fun` is praznik::MIM or
# praznik::MRMR. Returns the common contract: a length-p score (the per-feature MI
# criterion, -Inf for features the greedy search never ranked) and `selected` (the
# top-prop indices, or NULL as a pure ranker when prop is absent).
mi_select <- function(fun, who, time, status, X, prop, q) {
  if (!requireNamespace("praznik", quietly = TRUE)) {
    stop(sprintf("%s requires the 'praznik' package; install it with install.packages('praznik').", who),
         call. = FALSE)
  }
  s <- mi_substrate(time, status, X, q = q)
  k <- if (is.null(prop)) s$p else {
    stopifnot("`prop` must be a single number in (0, 1]" =
                length(prop) == 1L && is.numeric(prop) && prop > 0 && prop <= 1)
    max(1L, round(prop * s$p))
  }
  r <- fun(s$Xf, s$yf, k = k)
  score <- rep(-Inf, s$p)
  score[as.integer(r$selection)] <- as.numeric(r$score)
  selected <- if (is.null(prop)) NULL else as.integer(r$selection)
  list(score = score, selected = selected)
}

#' Survival mutual-information maximisation (direct MI baseline)
#'
#' Ranks features by their mutual information with survfcbf's own pseudo-target:
#' the martingale residual of a null Cox model, discretized into
#' equal-frequency bins exactly as [survfcbf_select()] does, then `praznik::MIM`.
#' Sharing the exact target transform and discretization isolates the scientific
#' question, since MIM and FCBF then differ only in the selection rule (plain
#' relevance ranking versus approximate-Markov-blanket predominance).
#'
#' @details
#' The pseudo-target and the per-feature cutpoints are fit on the supplied
#' (training) fold only, so the ranking carries no leakage. MIM is pure relevance,
#' so the score is order-independent mutual information; with `prop = NULL` a full
#' ranking is returned (`selected = NULL`), otherwise the top-`prop` features.
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. `NULL`
#'   (default) returns a pure ranking with `selected = NULL`.
#' @param q Discretization bins. `NULL` (default) uses Bommert's rule, matching
#'   [survfcbf_select()].
#' @param seed Unused (MIM is deterministic); present for interface uniformity.
#'
#' @return A list with `score` (per-feature mutual information, length `p`) and
#'   `selected` (top-`prop` indices, or `NULL`).
#'
#' @seealso [mrmr_surv_select()], [survfcbf_select()], [cut_ranking()]
#' @export
mim_surv_select <- function(time, status, X, prop = NULL, q = NULL, seed = NULL) {
  mi_select(praznik::MIM, "mim_surv_select", time, status, X, prop, q)
}

#' Survival mRMR (direct MI baseline)
#'
#' Ranks features by minimum-redundancy maximum-relevance on survfcbf's own
#' pseudo-target and discretization, via `praznik::MRMR`. Like
#' [mim_surv_select()] it shares the exact martingale target and equal-frequency
#' cutpoints with [survfcbf_select()], so the contrast against FCBF is purely the
#' selection rule: mRMR averages pairwise redundancy, FCBF prunes by dominance with
#' a Markov-blanket guarantee. mRMR was the family member that lost to a variance
#' filter in Bommert et al. (2022).
#'
#' @details
#' Fit on the training fold only. The greedy mRMR search is `O(p * k)`, so on wide
#' data the benchmark should pass a `prop` rather than rank all `p`. The score is
#' the mRMR criterion at each feature's selection step (`-Inf` for features never
#' reached by the search).
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. `NULL`
#'   (default) ranks all `p` (`selected = NULL`).
#' @param q Discretization bins. `NULL` (default) uses Bommert's rule.
#' @param seed Unused (mRMR is deterministic); present for interface uniformity.
#'
#' @return A list with `score` (per-feature mRMR criterion, length `p`) and
#'   `selected` (top-`prop` indices, or `NULL`).
#'
#' @seealso [mim_surv_select()], [survfcbf_select()], [cut_ranking()]
#' @export
mrmr_surv_select <- function(time, status, X, prop = NULL, q = NULL, seed = NULL) {
  mi_select(praznik::MRMR, "mrmr_surv_select", time, status, X, prop, q)
}
