#' Penalised Cox feature selection via glmnet (embedded baseline)
#'
#' Fits a penalised Cox model with `glmnet::cv.glmnet(family = "cox")`, the
#' embedded L1/L2 comparator. The penalty is `alpha = 1` for the
#' lasso (a sparse, self-selecting set) and `alpha = 0` for ridge (a dense ranker);
#' the regularisation `lambda` is chosen by inner cross-validation on the training
#' fold, so nothing leaks from held-out data.
#'
#' @details
#' Predictors are standardised internally by glmnet on the training matrix. For a
#' Cox model there is no intercept, so the coefficient vector maps one-to-one to
#' the `p` features. `score` is the absolute coefficient at the chosen `lambda`.
#' For `alpha > 0` the self-selected set is the nonzero coefficients; ridge
#' (`alpha = 0`) has no intrinsic sparsity and is a pure ranker (`selected = NULL`)
#' unless `prop` is given. `cv.glmnet` assigns CV folds at random, so `seed` makes
#' the call reproducible.
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. `NULL`
#'   (default) returns the nonzero set (lasso/elastic net) or `NULL` (ridge).
#' @param alpha Elastic-net mixing, `1` = lasso (default), `0` = ridge.
#' @param s Which cross-validated `lambda` to use, `"lambda.min"` (default) or
#'   `"lambda.1se"`.
#' @param nfolds Inner CV folds for `cv.glmnet`. Default `10`.
#' @param seed Integer or `NULL`. Seed for the reproducible inner CV.
#'
#' @return A list with `score` (absolute Cox coefficients, length `p`), `selected`
#'   (nonzero set, top-`prop`, or `NULL`), and `coef` (the signed coefficients).
#'
#' @seealso [rsf_select()], [cut_ranking()]
#' @export
cox_glmnet_select <- function(time, status, X, prop = NULL, alpha = 1,
                              s = c("lambda.min", "lambda.1se"), nfolds = 10, seed = NULL) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("cox_glmnet_select requires the 'glmnet' package; install it with install.packages('glmnet').",
         call. = FALSE)
  }
  s <- match.arg(s)
  stopifnot(
    "`X` must be a numeric matrix"                 = is.matrix(X) && is.numeric(X),
    "`X`, `time`, `status` must share nrow/length" = nrow(X) == length(time) && length(time) == length(status),
    "`alpha` must be a single number in [0, 1]"    = length(alpha) == 1L && is.numeric(alpha) && alpha >= 0 && alpha <= 1
  )
  if (is.logical(status)) status <- as.integer(status)
  y <- survival::Surv(time, status)
  cv <- with_seed(seed, glmnet::cv.glmnet(X, y, family = "cox", alpha = alpha,
                                          nfolds = nfolds, standardize = TRUE))
  co <- as.numeric(stats::coef(cv, s = s))           # length p, no intercept for a Cox model
  score <- abs(co)
  selected <- if (!is.null(prop)) {
    cut_ranking(score, prop = prop)
  } else if (alpha > 0) {
    which(co != 0)                                    # lasso / elastic net: nonzero = sparse self-selection
  } else {
    NULL                                              # ridge: no intrinsic sparsity, a pure ranker
  }
  list(score = score, selected = selected, coef = co)
}

#' Random survival forest feature selection via ranger (embedded baseline)
#'
#' Ranks features by random-survival-forest permutation importance using
#' `ranger`, a non-linear embedded comparator. It is a pure
#' ranker: the importance is the `score` and the top-`prop` features are
#' `selected`.
#'
#' @details
#' The forest is fit on the training fold only and the importance is read from it,
#' so there is no leakage. `ranger` is seeded through its own `seed` argument for
#' reproducibility. Permutation importance is the default as it is less biased
#' towards high-cardinality features than impurity importance.
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param prop Optional proportion of features to keep, in `(0, 1]`. `NULL`
#'   (default) returns a pure ranking (`selected = NULL`).
#' @param num.trees Number of trees. Default `500`.
#' @param importance Importance mode passed to `ranger`. Default `"permutation"`.
#' @param seed Integer or `NULL`. Seed passed to `ranger` for reproducibility.
#'
#' @return A list with `score` (per-feature importance, length `p`) and `selected`
#'   (top-`prop` indices, or `NULL`).
#'
#' @seealso [cox_glmnet_select()], [cut_ranking()]
#' @export
rsf_select <- function(time, status, X, prop = NULL, num.trees = 500,
                       importance = "permutation", seed = NULL) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("rsf_select requires the 'ranger' package; install it with install.packages('ranger').",
         call. = FALSE)
  }
  stopifnot(
    "`X` must be a numeric matrix"                 = is.matrix(X) && is.numeric(X),
    "`X`, `time`, `status` must share nrow/length" = nrow(X) == length(time) && length(time) == length(status)
  )
  if (is.logical(status)) status <- as.integer(status)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))   # ranger needs covariate names
  y <- survival::Surv(time, status)
  rf <- ranger::ranger(x = X, y = y, num.trees = num.trees, importance = importance,
                       seed = if (is.null(seed)) NULL else as.integer(seed), verbose = FALSE)
  score <- as.numeric(rf$variable.importance)
  selected <- if (is.null(prop)) NULL else cut_ranking(score, prop = prop)
  list(score = score, selected = selected)
}
