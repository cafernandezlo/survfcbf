# internal: correlation-block grouping by union-find on |cor(X)| > rho_t. Features
# whose absolute pairwise correlation exceeds the threshold share a group; the rest
# are singletons. Fit on the training fold (the correlation matrix is learned here).
correlation_blocks <- function(X, rho_t) {
  p <- ncol(X)
  parent <- seq_len(p)
  find <- function(i) {
    root <- i
    while (parent[root] != root) root <- parent[root]
    while (parent[i] != root) { nxt <- parent[i]; parent[i] <<- root; i <- nxt }   # path compression
    root
  }
  C <- abs(stats::cor(X))
  for (i in seq_len(p - 1L)) {
    hits <- which(C[i, (i + 1L):p] > rho_t)          # columns > i correlated with i (NA dropped)
    for (jj in hits) {
      a <- find(i); b <- find(i + jj)
      if (a != b) parent[a] <- b                     # local update; find() does the <<- path compression
    }
  }
  as.integer(factor(vapply(seq_len(p), find, integer(1))))   # contiguous group ids 1..G
}

# internal: V_b per group, the largest lambda at which the group is nonzero (its
# entry point on the path); 0 if the group never enters. beta is features x lambda.
group_entry_lambda <- function(beta, group, lambda) {
  vapply(sort(unique(group)), function(g) {
    rows <- which(group == g)
    nz <- apply(abs(beta[rows, , drop = FALSE]) > 0, 2L, any)
    if (any(nz)) lambda[which(nz)[1L]] else 0
  }, numeric(1))
}

# internal: permutation-assisted group-Lasso core (Le et al. 2024 stage 2), run on a
# given feature matrix. Groups features by correlation, then over K rounds fits a
# group-Lasso Cox model alongside a row-permuted pseudo mirror and keeps groups whose
# path entry lambda V_b beats the strongest pseudo group in more than tau of rounds.
# Assumes the caller has set the RNG stream (it draws permutations); deterministic
# grpsurv otherwise. Returns per-feature score (group frequency) and selected indices.
pagl_core <- function(time, status, X, rho_t, K, tau, penalty) {
  n <- nrow(X)
  p <- ncol(X)
  groups <- correlation_blocks(X, rho_t)
  G <- max(groups)
  group_aug <- c(groups, groups + G)                 # real groups, then their pseudo mirrors
  y <- survival::Surv(time, status)

  acc <- integer(G)
  for (k in seq_len(K)) {
    Xaug <- cbind(X, X[sample.int(n), , drop = FALSE])   # pseudo = row-permuted features (preserves correlation)
    colnames(Xaug) <- as.character(seq_len(2L * p))      # unique names so grpsurv accepts them
    fit <- grpreg::grpsurv(Xaug, y, group = group_aug, penalty = penalty)
    Vb <- group_entry_lambda(fit$beta, group_aug, fit$lambda)
    acc <- acc + as.integer(Vb[seq_len(G)] > max(Vb[(G + 1L):(2L * G)]))
  }
  group_frequency <- acc / K
  list(group = groups, group_frequency = group_frequency,
       score = group_frequency[groups],
       selected = which(groups %in% which(group_frequency > tau)))
}

# internal: a single base selector (stage 1 ensemble member) -> selected column
# indices (a set). Rankers (cox, mim, mrmr) keep up to n_select features; the lasso
# self-selects its nonzero set. All fit on the supplied (inner-fold) data only.
base_selector_set <- function(which_sel, time, status, X, n_select, cox_alpha) {
  p <- ncol(X)
  prop <- min(1, n_select / p)
  out <- switch(which_sel,
    lasso = tryCatch({ r <- cox_glmnet_select(time, status, X, alpha = 1); if (is.null(r$selected)) integer(0) else r$selected },
                     error = function(e) integer(0)),
    cox   = univariate_cox_set(time, status, X, k = n_select, alpha = cox_alpha),
    mim   = tryCatch(mim_surv_select(time, status, X, prop = prop)$selected,  error = function(e) integer(0)),
    mrmr  = tryCatch(mrmr_surv_select(time, status, X, prop = prop)$selected, error = function(e) integer(0))
  )
  as.integer(out)
}

# internal: univariate Cox screening (Le et al. "Cox" selector) -- up to the top
# `k` features with Wald p-value below `alpha`, fit one feature at a time.
univariate_cox_set <- function(time, status, X, k = 25L, alpha = 0.05) {
  y <- survival::Surv(time, status)
  pvals <- vapply(seq_len(ncol(X)), function(j) {
    if (stats::sd(X[, j]) == 0) return(1)
    fit <- tryCatch(survival::coxph(y ~ X[, j]), error = function(e) NULL)
    if (is.null(fit)) return(1)
    co <- summary(fit)$coefficients
    if (nrow(co) == 0L) 1 else co[1L, "Pr(>|z|)"]
  }, numeric(1))
  cand <- which(pvals < alpha)
  if (length(cand) == 0L) return(integer(0))
  cand[order(pvals[cand])][seq_len(min(k, length(cand)))]
}

# internal: stage 1 weighted-relative-frequency ensemble (He et al. 2022 scheme,
# Le et al. 2024 base selectors). For each selector and each of cv_folds folds, run
# it on the all-but-fold data and add 1/|selected| to each picked feature's WRF;
# convert each selector's WRF to ranks (higher WRF ranks higher), average ranks
# across selectors, and return the candidate set (any feature ever selected).
wrf_stage1 <- function(time, status, X, cv_folds, n_select, ensemble, cox_alpha) {
  n <- nrow(X)
  p <- ncol(X)
  selectors <- c("lasso", "cox", "mim", if (ensemble >= 2L) "mrmr")
  folds <- sample(rep_len(seq_len(cv_folds), n))
  wrf <- matrix(0, p, length(selectors))
  for (si in seq_along(selectors)) {
    for (f in seq_len(cv_folds)) {
      idx <- which(folds != f)                       # all-but-fold-f training data
      sel <- base_selector_set(selectors[si], time[idx], status[idx], X[idx, , drop = FALSE],
                               n_select, cox_alpha)
      if (length(sel) > 0L) wrf[sel, si] <- wrf[sel, si] + 1 / length(sel)
    }
  }
  ranks <- apply(wrf, 2L, function(w) rank(-w, ties.method = "average"))   # higher WRF -> lower rank
  list(candidate = which(rowSums(wrf) > 0),          # union of every selection
       mean_rank = rowMeans(ranks),
       wrf = rowSums(wrf))
}

#' Permutation-assisted group-Lasso ensemble (Le et al. 2024 baseline)
#'
#' The most elaborate baseline: the two-stage ensemble of Le et
#' al. (2024). Stage 1 aggregates several base selectors over a cross-validation
#' resampling into a candidate feature set by weighted relative frequency; stage 2
#' fits a group-Lasso Cox model on that candidate set with permutation-assisted
#' tuning, keeping groups that consistently beat permuted pseudo-feature controls.
#'
#' @details
#' Stage 1 follows the He et al. (2022) weighted-relative-frequency scheme that Le
#' et al. adapt: for each base selector and each of `cv_folds` folds, the selector
#' runs on the all-but-fold data and each picked feature accrues `1 / (number
#' selected)`. The base selectors are the survival-adapted members of Le et al.'s
#' Ensemble 2: the Cox lasso ([cox_glmnet_select()]), univariate Cox screening (up
#' to `n_select` features with Wald `p < cox_alpha`), [mim_surv_select()] and (when
#' `ensemble = 2`) [mrmr_surv_select()]. The candidate set is every feature selected
#' at least once.
#'
#' Stage 2 (the permutation-assisted group lasso) groups the candidates into
#' correlation blocks (`|cor| > rho_t`, default `0.75`), and over `K` rounds fits a
#' group-Lasso Cox model (`grpreg::grpsurv`) on the candidates plus a row-permuted
#' pseudo mirror (a knockoff-style null that preserves the feature correlation
#' structure but breaks the outcome link). A real group is selected in a round when
#' its path entry `lambda` (`V_b`, the largest `lambda` at which it is nonzero)
#' beats the strongest pseudo group; groups kept in more than `tau` of rounds form
#' the set. Everything is fit on the training fold; `seed` makes it reproducible.
#'
#' Choices the papers leave open, set here as documented defaults: cross-selector
#' aggregation by mean rank, the candidate set as the union of selections, the
#' pseudo control as a whole-matrix row permutation, and `V_b` as the path entry.
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. `NULL`
#'   (default) returns the features of the retained groups (the self-selected set).
#' @param rho_t Absolute-correlation threshold for the stage-2 grouping. Default
#'   `0.75`.
#' @param K Number of permutation rounds in stage 2. Default `50`.
#' @param tau Selection-frequency threshold over the rounds, in `(0, 1)`. Default
#'   `0.5`.
#' @param penalty Group penalty passed to `grpreg::grpsurv`. Default `"grLasso"`.
#' @param cv_folds Stage-1 cross-validation folds for the WRF. Default `5`.
#' @param n_select Features kept per ranking base selector (and the cap for the Cox
#'   screen) in stage 1. Default `25`.
#' @param ensemble `2` (default) uses Lasso + Cox + MIM + MRMR, `1` drops MRMR.
#' @param cox_alpha Wald p-value cutoff for the univariate Cox base selector.
#'   Default `0.05`.
#' @param seed Integer or `NULL`. Seed for the reproducible resampling and rounds.
#'
#' @return A list with `score` (per-feature stage-2 group selection frequency,
#'   length `p`, `0` for non-candidates), `selected` (retained-group features, or
#'   the top-`prop`), `candidate` (stage-1 candidate indices), `group` (length-`p`
#'   group id among candidates, `0` for non-candidates), and `group_frequency`.
#'
#' @seealso [carss_select()], [cox_glmnet_select()], [cut_ranking()]
#' @export
grplasso_ensemble_select <- function(time, status, X, prop = NULL,
                                     rho_t = 0.75, K = 50L, tau = 0.5, penalty = "grLasso",
                                     cv_folds = 5L, n_select = 25L, ensemble = 2L, cox_alpha = 0.05,
                                     seed = NULL) {
  need <- c("grpreg", "glmnet", "praznik")             # stage 2 (grpreg) + stage-1 base selectors
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0L) {
    stop(sprintf("grplasso_ensemble_select needs %s; install with install.packages(c(%s)).",
                 paste(miss, collapse = ", "), paste(sprintf("'%s'", miss), collapse = ", ")),
         call. = FALSE)
  }
  stopifnot(
    "`X` must be a numeric matrix"                 = is.matrix(X) && is.numeric(X),
    "`X`, `time`, `status` must share nrow/length" = nrow(X) == length(time) && length(time) == length(status),
    "`rho_t` must be a single number in (0, 1)"    = length(rho_t) == 1L && is.numeric(rho_t) && rho_t > 0 && rho_t < 1,
    "`K` must be a positive integer"               = length(K) == 1L && is.numeric(K) && is.finite(K) && K >= 1 && K == as.integer(K),
    "`tau` must be a single number in (0, 1)"      = length(tau) == 1L && is.numeric(tau) && tau > 0 && tau < 1,
    "`cv_folds` must be an integer >= 2"           = length(cv_folds) == 1L && is.numeric(cv_folds) && cv_folds >= 2 && cv_folds == as.integer(cv_folds),
    "`ensemble` must be 1 or 2"                    = length(ensemble) == 1L && ensemble %in% c(1L, 2L)
  )
  if (is.logical(status)) status <- as.integer(status)
  K <- as.integer(K)
  p <- ncol(X)

  res <- with_seed(seed, {
    s1 <- wrf_stage1(time, status, X, as.integer(cv_folds), n_select, as.integer(ensemble), cox_alpha)
    cand <- s1$candidate
    if (length(cand) < 2L) {
      list(candidate = cand, group = integer(0), group_frequency = numeric(0),
           score_local = rep(0, length(cand)), selected_local = seq_along(cand))
    } else {
      pg <- pagl_core(time, status, X[, cand, drop = FALSE], rho_t, K, tau, penalty)
      list(candidate = cand, group = pg$group, group_frequency = pg$group_frequency,
           score_local = pg$score, selected_local = pg$selected)
    }
  })

  cand <- res$candidate
  score <- numeric(p)
  group <- integer(p)
  if (length(cand) > 0L) {
    score[cand] <- res$score_local
    if (length(res$group) > 0L) group[cand] <- res$group
  }
  selected <- if (!is.null(prop)) cut_ranking(score, prop = prop) else cand[res$selected_local]
  list(score = score, selected = as.integer(selected), candidate = cand,
       group = group, group_frequency = res$group_frequency)
}
