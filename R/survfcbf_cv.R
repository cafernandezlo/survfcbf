#' survfcbf with inner-CV-tuned hyperparameters
#'
#' Tunes survfcbf's discretization `q` (and optionally the relevance threshold
#' `delta`) by inner cross-validation, choosing the configuration that optimises a
#' downstream survival learner's predictive performance on held-out inner folds.
#' This is the headline adaptive survfcbf configuration and, with a `delta_grid`,
#' the prediction-optimal "at-threshold" variant (the cardinality knob that keeps
#' the predominance pruning but relaxes the conservative permutation threshold so a
#' Cox learner has enough features). Standard nested-CV tuning, so the selection is
#' learner-coupled by construction.
#'
#' @details
#' Inner CV is stratified on the event indicator and seeded. For each candidate
#' (`q`, `delta`) it runs [survfcbf_select()] on each inner-train fold, fits a
#' clone of `learner` on the selected features, and scores the inner-test fold with
#' `inner_measure` (Uno's C by default). Candidates whose selection is too small
#' for the learner score `NA` and are skipped, so the tuner naturally avoids
#' degenerate configurations. The best candidate is refit on the full training
#' data. With `delta_grid = NULL` only `q` is tuned and `delta` is the permutation
#' null per candidate (headline); a `delta_grid` tunes the threshold too (lower
#' `delta` keeps more predominant features).
#'
#' Fit on the training fold only and seeded. Requires `mlr3proba` and a survival
#' learner.
#'
#' @param time,status,X Training-fold outcome and feature matrix.
#' @param learner An `mlr3` survival learner used to score candidates (and, by the
#'   caller, downstream).
#' @param q_grid Candidate discretization bin counts. Default `c(4, 7, 10)`.
#' @param delta_grid Candidate relevance thresholds, or `NULL` (default) to take
#'   `delta` from the permutation null per candidate (tune `q` only).
#' @param inner_folds Inner CV folds. Default `5`.
#' @param B Permutations for the permutation-null `delta`. Default `200`.
#' @param target_method,relevance Passed to [survfcbf_select()].
#' @param inner_measure An `mlr3` survival measure to optimise. `NULL` uses Uno's C.
#' @param seed Integer or `NULL`. Seed for the inner CV and the refit.
#'
#' @return A list with `selected` (predominant indices at the chosen config), the
#'   chosen `q` and `delta`, `inner_score` (best mean inner-CV score), `relevance`,
#'   and `cv` (the per-candidate score table).
#'
#' @seealso [survfcbf_select()], [nested_cv_select()]
#' @export
survfcbf_cv_select <- function(time, status, X, learner, q_grid = c(4, 7, 10),
                               delta_grid = NULL, inner_folds = 5L, B = 200L,
                               target_method = "martingale", relevance = "su",
                               inner_measure = NULL, seed = NULL) {
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("survfcbf_cv_select requires the 'mlr3proba' package.", call. = FALSE)
  }
  stopifnot(
    "`X` must be a numeric matrix"                 = is.matrix(X) && is.numeric(X),
    "`X`, `time`, `status` must share nrow/length" = nrow(X) == length(time) && length(time) == length(status),
    "`inner_folds` must be an integer >= 2"        = length(inner_folds) == 1L && inner_folds >= 2 && inner_folds == as.integer(inner_folds)
  )
  if (is.null(inner_measure)) inner_measure <- mlr3::msr("surv.cindex", weight_meth = "G2")
  minimize <- isTRUE(inner_measure$minimize)
  cand <- expand.grid(q = q_grid, delta = if (is.null(delta_grid)) NA_real_ else delta_grid)

  with_seed(seed, {
    folds <- stratified_folds(status, as.integer(inner_folds))
    score_candidate <- function(q, d) {
      fs <- vapply(seq_len(inner_folds), function(f) {
        tr <- which(folds != f)
        te <- which(folds == f)
        sel <- tryCatch(
          survfcbf_select(time[tr], status[tr], X[tr, , drop = FALSE], q = q, B = B,
                          target_method = target_method, relevance = relevance,
                          delta = if (is.na(d)) NULL else d)$selected,
          error = function(e) integer(0))
        if (length(sel) < 1L) return(NA_real_)               # too few features for the learner
        task <- surv_task(time, status, X[, sel, drop = FALSE])
        rid <- task$row_ids
        lf <- learner$clone(deep = TRUE)
        pred <- tryCatch({ lf$train(task, row_ids = rid[tr]); lf$predict(task, row_ids = rid[te]) },
                         error = function(e) NULL)
        if (is.null(pred)) NA_real_ else pred$score(inner_measure, task = task, train_set = rid[tr])
      }, numeric(1))
      mean(fs, na.rm = TRUE)
    }
    scores <- vapply(seq_len(nrow(cand)), function(i) score_candidate(cand$q[i], cand$delta[i]), numeric(1))
    best <- if (all(is.na(scores))) 1L else if (minimize) which.min(scores) else which.max(scores)
    final <- survfcbf_select(time, status, X, q = cand$q[best], B = B,
                             target_method = target_method, relevance = relevance,
                             delta = if (is.na(cand$delta[best])) NULL else cand$delta[best])
    list(selected = final$selected, q = cand$q[best], delta = final$delta,
         inner_score = scores[best], relevance = final$relevance,
         cv = data.frame(cand, score = scores, row.names = NULL))
  })
}
