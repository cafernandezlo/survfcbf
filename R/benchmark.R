#' Build an mlr3proba survival task from the matrix form
#'
#' Convenience wrapper turning the `(time, status, X)` form used throughout the
#' package into an `mlr3proba::TaskSurv` for the benchmark harness.
#' Feature names are taken from `colnames(X)` or generated.
#'
#' @param time,status Right-censored outcome.
#' @param X Numeric `n` by `p` feature matrix.
#' @param id Task id. Default `"surv"`.
#'
#' @return An `mlr3proba::TaskSurv`.
#' @seealso [nested_cv_select()]
#'
#' @examples
#' if (requireNamespace("mlr3proba", quietly = TRUE)) {
#'   sim <- simulate_dataset(default_params(
#'     p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'     n_distractor_subblocks = 2L, distractor_size = 4L,
#'     calib_n_presample = 500, seed = 1))
#'   task <- surv_task(sim$time, sim$status, sim$X)
#'   task$nrow
#'   length(task$feature_names)
#' }
#' @export
surv_task <- function(time, status, X, id = "surv") {
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("surv_task requires the 'mlr3proba' package.", call. = FALSE)
  }
  stopifnot("`X` must be a numeric matrix" = is.matrix(X) && is.numeric(X),
            "`X`, `time`, `status` must share nrow/length" =
              nrow(X) == length(time) && length(time) == length(status))
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  d <- data.frame(time = as.numeric(time), status = as.integer(status), X, check.names = FALSE)
  mlr3proba::as_task_surv(d, time = "time", event = "status", id = id)
}

# internal: fold assignment 1..k stratified on the event indicator, so each fold
# carries a similar event rate. Seeded for reproducibility.
stratified_folds <- function(status, k, seed = NULL) {
  with_seed(seed, {
    fold <- integer(length(status))
    for (lvl in unique(status)) {
      idx <- which(status == lvl)
      fold[idx] <- sample(rep_len(seq_len(k), length(idx)))
    }
    fold
  })
}

# internal: pull the (time, status, X, feature names) back out of a TaskSurv. Uses
# the Surv object for the outcome, whose structure ([, 1] time, [, 2] status) is
# fixed, so it is robust to the task's target column names (real-data tasks in 7.6).
extract_surv <- function(task) {
  fn <- task$feature_names
  surv <- task$truth()                          # Surv for all rows, in row order
  list(time = as.numeric(surv[, 1L]), status = as.integer(surv[, 2L]),
       X = as.matrix(task$data(cols = fn)), features = fn)
}

#' Leakage-free nested cross-validation of a selector and learner
#'
#' The benchmark engine. Over `outer_folds` folds
#' (stratified on the event indicator, seeded), each outer-train fold has its
#' features selected by `select_fn`, the task is reduced to the selected features,
#' the downstream `learner` is trained on that reduced training fold, and a
#' prediction is produced on the held-out outer-test fold. Selection and the learner
#' see only training-fold rows, so there is no leakage; the returned predictions and
#' selected sets feed the metric, stability and parsimony layers.
#'
#' @details
#' `select_fn` is a function `(time, status, X)` returning the selected column
#' indices, so any selector (survfcbf or a baseline) is plugged in by a thin
#' closure that bakes in its hyperparameters. Each fold runs inside its own seeded
#' RNG stream, so the selector and the learner are jointly reproducible; the
#' selectors and the learner are therefore called with their own seeds left unset.
#' A fold whose selector returns nothing is recorded as an empty fold (cardinality
#' `0`, `prediction` `NULL`, `error` `"empty selection"`) rather than silently kept
#' as all features, and a learner that fails on a degenerate selection is caught and
#' recorded in `error` rather than crashing the benchmark.
#'
#' @param task An `mlr3proba::TaskSurv` (see [surv_task()]).
#' @param select_fn A function `(time, status, X)` returning integer column
#'   indices, e.g. `function(ti, st, X) carss_select(ti, st, X)$selected`.
#' @param learner An `mlr3` survival learner, e.g.
#'   `mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")`. Prefer
#'   `s = "lambda.min"` for the cross-validated ridge Cox: under the learner's
#'   default `"lambda.1se"` the one-standard-error rule picks the first lambda of the
#'   ridge path whenever events are scarce relative to `p`, the fitted model is then
#'   the null model, and its predicted distribution is the same for every panel (see
#'   [score_benchmark()], `null_fit`).
#' @param outer_folds Number of outer CV folds. Default `10`.
#' @param seed Integer or `NULL`. Master seed for the folds and per-fold streams.
#' @param rng How per-fold seeds are derived from `seed`. `"sequential"` (default)
#'   derives arithmetic per-fold seeds, the historical behaviour, kept as the
#'   default so results produced with it reproduce exactly. `"stream"` spaces one
#'   independent L'Ecuyer-CMRG stream per fold with [parallel::nextRNGStream()],
#'   the same construction the simulator uses ([run_replicates()]) and the
#'   recommended choice for new studies, since sequential integer seeds do not
#'   guarantee independent random number streams across folds. Requires a
#'   non-`NULL` `seed`. The fold assignment itself is identical under both, so
#'   only the selector's and learner's random draws differ.
#'
#' @return An object of class `survfcbf_benchmark`: a list with `folds` (one record
#'   per outer fold: `prediction` (`NULL` if the learner failed), `train_set`,
#'   `test_set`, `selected`, `features`, `cardinality`, `select_time`, `error`),
#'   `fold_id`, the `task`, and `outer_folds`.
#'
#' @seealso [surv_task()]
#'
#' @examples
#' if (requireNamespace("mlr3proba", quietly = TRUE)) {
#'   sim <- simulate_dataset(default_params(
#'     p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'     n_distractor_subblocks = 2L, distractor_size = 4L,
#'     calib_n_presample = 500, seed = 1))
#'   task <- surv_task(sim$time, sim$status, sim$X)
#'
#'   # any (time, status, X) -> column indices function is a selector
#'   top5 <- function(time, status, X)
#'     survfcbf_select(time, status, X, q = 2, top = 5)$selected
#'
#'   bm <- nested_cv_select(task, top5, mlr3::lrn("surv.coxph"),
#'                          outer_folds = 3, seed = 1)
#'   score_benchmark(bm)$summary   # Uno's C (and truncated), Harrell's C, IBS
#' }
#' @export
nested_cv_select <- function(task, select_fn, learner, outer_folds = 10L, seed = NULL,
                             rng = c("sequential", "stream")) {
  rng <- match.arg(rng)
  if (rng == "stream" && is.null(seed)) {
    stop("`rng = \"stream\"` needs a `seed` to derive the per-fold streams from.", call. = FALSE)
  }
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("nested_cv_select requires the 'mlr3proba' package.", call. = FALSE)
  }
  stopifnot(
    "`select_fn` must be a function"               = is.function(select_fn),
    "`outer_folds` must be an integer >= 2"        = length(outer_folds) == 1L && outer_folds >= 2 && outer_folds == as.integer(outer_folds)
  )
  outer_folds <- as.integer(outer_folds)
  d <- extract_surv(task)
  rids <- task$row_ids
  folds <- stratified_folds(d$status, outer_folds, seed)
  streams <- if (rng == "stream") fold_rng_streams(seed, outer_folds) else NULL

  res <- lapply(seq_len(outer_folds), function(f) {
    tr <- which(folds != f)
    te <- which(folds == f)
    ## "sequential" keeps the historical arithmetic per-fold seeds so existing results
    ## reproduce bit for bit; "stream" evaluates the fold under its own L'Ecuyer-CMRG
    ## stream state instead (see fold_rng_streams), which is the package invariant.
    fold_seed <- if (is.null(seed) || rng == "stream") NULL else
      as.integer((as.numeric(seed) * 1000 + f) %% 2147483647)
    with_rng_state(if (rng == "stream") streams[[f]] else NULL, with_seed(fold_seed, {
      t0 <- proc.time()[["elapsed"]]
      sel <- as.integer(select_fn(d$time[tr], d$status[tr], d$X[tr, , drop = FALSE]))
      select_time <- proc.time()[["elapsed"]] - t0
      base <- list(train_set = rids[tr], test_set = rids[te], select_time = select_time)
      if (length(sel) == 0L) {
        # the selector chose nothing: record an honest empty fold (no model), do not
        # mask it as a no-filter selection (that would inflate parsimony and predict)
        c(base, list(prediction = NULL, selected = integer(0), features = character(0),
                     cardinality = 0L, error = "empty selection"))
      } else {
        feats <- d$features[sel]
        tsel <- task$clone()$select(feats)
        lrn_f <- learner$clone(deep = TRUE)
        err <- NULL
        pred <- tryCatch({                                # a degenerate selection can still break a learner
          lrn_f$train(tsel, row_ids = rids[tr])           # (e.g. glmnet needs >= 2 features); record the
          lrn_f$predict(tsel, row_ids = rids[te])         # failure rather than crash the whole benchmark
        }, error = function(e) { err <<- conditionMessage(e); NULL })
        c(base, list(prediction = pred, selected = sel, features = feats,
                     cardinality = length(sel), error = err))
      }
    }))
  })

  structure(list(folds = res, fold_id = folds, task = task, outer_folds = outer_folds),
            class = "survfcbf_benchmark")
}

#' Default survival measures for the benchmark
#'
#' The discrimination and calibration measures scored on each fold:
#' Uno's C (IPCW-corrected, untruncated), Harrell's C (for continuity with the
#' wider literature), and the integrated Brier score (IPCW Graf, discrimination
#' plus calibration, for comparability with Bommert). The primary discrimination
#' metric, Uno's C truncated at a horizon learned on each training fold, is added
#' by [score_benchmark()] through its `truncation` argument, since its horizon
#' differs between folds and cannot live in a fixed measure.
#'
#' @return A named list of `mlr3` survival measures (`uno_c`, `harrell_c`, `ibs`).
#' @seealso [score_benchmark()]
#' @export
default_surv_measures <- function() {
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("default_surv_measures requires the 'mlr3proba' package.", call. = FALSE)
  }
  list(
    uno_c     = mlr3::msr("surv.cindex", weight_meth = "G2", id = "uno_c"),
    harrell_c = mlr3::msr("surv.cindex", id = "harrell_c"),
    ibs       = mlr3::msr("surv.graf", id = "ibs")
  )
}

#' Score a nested-CV benchmark with Monte Carlo standard errors
#'
#' Scores each outer fold of a [nested_cv_select()] result with survival measures
#' and aggregates across folds. The IPCW measures (Uno's C, the
#' integrated Brier) are given the fold's training rows so the censoring
#' distribution is estimated from training data, with no leakage. Folds whose
#' learner failed (`prediction` is `NULL`) are skipped and counted, not hidden.
#'
#' @details
#' A fold is flagged as a null fit (`null_fit`) when its predicted survival curves
#' are the same for every held-out row, i.e. the trained model ignored the selected
#' features. Such a fold is still scored, since a learner that settles on the null
#' model is a legitimate outcome, but its integrated Brier score is the score of the
#' training baseline and cannot separate panels. A cross-validated ridge Cox tuned
#' by the one-standard-error rule produces exactly this when the rule picks the first
#' lambda of the path: the coefficients are of order `1e-37`, so the predicted
#' distribution is the baseline for any panel while the C-index still ranks on the
#' vanishing linear predictor. A selector-independent IBS is the symptom to look
#' for; `n_null_fit` makes it visible.
#'
#' Untruncated, Uno's C weights each comparable pair by the inverse squared
#' censoring survival estimated on the training rows. A test event later than the
#' last censoring time of the training fold meets an estimate of zero, its pairs
#' take a capped weight of `1 / eps^2`, and a handful of pairs can carry nearly all
#' the weight of the fold. `truncation` therefore also scores Uno's C restricted to
#' times up to `t_max`, the `truncation` quantile of the observed times of the
#' fold's training rows (Uno et al. 2011), as `uno_c_trunc`. The horizon is learned
#' on the training rows only, and is reported per fold in `per_fold$t_max`.
#'
#' @param bm A `survfcbf_benchmark` from [nested_cv_select()].
#' @param measures A named list of `mlr3` survival measures. `NULL` (default) uses
#'   [default_surv_measures()].
#' @param truncation Quantile of the training fold's observed times at which the
#'   truncated Uno's C (`uno_c_trunc`) is cut, a number in `(0, 1]`. Default `0.8`.
#'   `NULL` skips it.
#'
#' @return A list with `per_fold` (a data.frame: `fold`, `cardinality`,
#'   `select_time`, `null_fit`, one column per measure, then `uno_c_trunc` and
#'   `t_max` unless `truncation` is `NULL`), `summary` (per measure:
#'   `mean`, `se` = `sd / sqrt(folds)`, `n_folds` scored), `n_failed` (folds with no
#'   prediction) and `n_null_fit` (scored folds whose predicted distribution does
#'   not vary across the held-out rows).
#'
#' @references Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the
#'   C-statistics for evaluating overall adequacy of risk prediction procedures with
#'   censored survival data. *Statistics in Medicine* 30(10):1105-1117.
#'
#' @seealso [nested_cv_select()], [default_surv_measures()]
#' @export
score_benchmark <- function(bm, measures = NULL, truncation = 0.8) {
  if (!requireNamespace("mlr3proba", quietly = TRUE)) {
    stop("score_benchmark requires the 'mlr3proba' package.", call. = FALSE)
  }
  stopifnot(
    "`bm` must be a survfcbf_benchmark" = inherits(bm, "survfcbf_benchmark"),
    "`truncation` must be NULL or a single number in (0, 1]" = is.null(truncation) ||
      (is.numeric(truncation) && length(truncation) == 1L && !is.na(truncation) &&
         truncation > 0 && truncation <= 1)
  )
  if (is.null(measures)) measures <- default_surv_measures()
  mids <- vapply(measures, function(m) m$id, character(1))
  if (!is.null(truncation) && "uno_c_trunc" %in% mids) {
    stop("a measure already has the id `uno_c_trunc`: rename it or set `truncation = NULL`.",
         call. = FALSE)
  }
  ## the truncation horizon of each fold, from its training rows only
  t_max <- if (is.null(truncation)) NULL else vapply(bm$folds, function(fd)
    unname(stats::quantile(bm$task$truth(rows = fd$train_set)[, 1L], truncation)), numeric(1))
  ids <- c(mids, if (!is.null(truncation)) "uno_c_trunc")

  scored <- lapply(seq_along(bm$folds), function(f) {
    fd <- bm$folds[[f]]
    if (is.null(fd$prediction)) return(stats::setNames(rep(NA_real_, length(ids)), ids))
    ms <- measures
    if (!is.null(truncation)) {
      ms <- c(ms, list(mlr3::msr("surv.cindex", weight_meth = "G2", t_max = t_max[f],
                                 id = "uno_c_trunc")))
    }
    v <- fd$prediction$score(ms, task = bm$task, train_set = fd$train_set)
    v[ids]                                           # fix order and names
  })

  per_fold <- data.frame(
    fold = seq_along(bm$folds),
    cardinality = vapply(bm$folds, function(fd) as.numeric(fd$cardinality), numeric(1)),
    select_time = vapply(bm$folds, function(fd) as.numeric(fd$select_time), numeric(1)),
    null_fit = vapply(bm$folds, function(fd)
      if (is.null(fd$prediction)) NA else constant_distr(fd$prediction), logical(1)),
    do.call(rbind, scored), row.names = NULL, check.names = FALSE
  )
  if (!is.null(truncation)) per_fold$t_max <- t_max

  summary <- do.call(rbind, lapply(ids, function(m) {
    x <- per_fold[[m]]
    x <- x[is.finite(x)]
    data.frame(metric = m, mean = if (length(x)) mean(x) else NA_real_,
               se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_,
               n_folds = length(x), row.names = NULL)
  }))

  list(per_fold = per_fold, summary = summary,
       n_failed = sum(vapply(bm$folds, function(fd) is.null(fd$prediction), logical(1))),
       n_null_fit = sum(per_fold$null_fit, na.rm = TRUE))
}

# internal: TRUE when a prediction's survival curves are identical across its rows
# (the model ignored its features), NA when there is no survival matrix to inspect
# (a crank-only learner, or a single held-out row).
constant_distr <- function(pred, tol = sqrt(.Machine$double.eps)) {
  S <- pred$data$distr
  if (!is.matrix(S) || nrow(S) < 2L) return(NA)
  all(abs(S - rep(S[1L, ], each = nrow(S))) <= tol)
}

#' Feature-similarity matrix for the stability metric
#'
#' Builds the `p` by `p` absolute-correlation matrix that
#' [benchmark_stability()] uses to credit a selector for choosing a correlated
#' substitute instead of the identical feature.
#'
#' @details
#' With `scope = "full"` the matrix is `abs(cor(X))` over every row of the task.
#' With `scope = "train"` it is the element-wise mean of the per-fold
#' training-set correlation matrices, computed from the benchmark's stored fold
#' assignment, so no held-out row contributes. The matrix is unsupervised either
#' way (no outcome enters it) and is only used to score agreement between
#' already-selected sets, never to select; `"train"` exists for studies where any
#' use of held-out rows must be ruled out. Correlations that are undefined on a
#' training subset (a feature constant within it) are set to `0`, and the
#' diagonal to `1`.
#'
#' @param bm A `survfcbf_benchmark` from [nested_cv_select()].
#' @param scope `"full"` (all rows) or `"train"` (mean over per-fold training
#'   rows).
#'
#' @return A `p` by `p` numeric matrix of absolute correlations.
#'
#' @seealso [benchmark_stability()], [benchmark_selectors()]
#' @export
similarity_matrix <- function(bm, scope = c("full", "train")) {
  scope <- match.arg(scope)
  stopifnot("`bm` must be a survfcbf_benchmark" = inherits(bm, "survfcbf_benchmark"))
  X <- extract_surv(bm$task)$X
  m <- if (scope == "full") {
    abs(stats::cor(X))
  } else {
    stopifnot("`scope = \"train\"` needs the benchmark's fold assignment (`bm$fold_id`)" =
                !is.null(bm$fold_id) && length(bm$fold_id) == nrow(X))
    acc <- NULL                                            # running sum: one p x p accumulator, so the peak
    ks <- sort(unique(bm$fold_id))                         # footprint is ~2 matrices, not folds-many
    for (k in ks) {
      mk <- abs(stats::cor(X[bm$fold_id != k, , drop = FALSE]))
      mk[is.na(mk)] <- 0                                   # feature constant within this training subset
      acc <- if (is.null(acc)) mk else acc + mk
    }
    acc / length(ks)
  }
  m[is.na(m)] <- 0
  diag(m) <- 1
  m
}

#' Selection stability, parsimony and runtime of a benchmark
#'
#' Summarises the selection quality of a [nested_cv_select()] result:
#' correlation-aware stability across the outer-fold selected sets, parsimony
#' (selected cardinality) and selection runtime.
#'
#' @details
#' Stability is Bommert's correlation-adjusted count-based measure (SMA-Count,
#' `stabm::stabilityIntersectionCount`) over the per-fold selected sets, with a
#' similarity matrix of absolute feature correlations thresholded at `theta`
#' (default `0.9`), so swapping a feature for a strongly correlated one is not
#' penalised. The chance correction uses sampling, so `seed` makes it reproducible.
#' Folds that selected nothing do not contribute to stability or parsimony (and are
#' counted in `n_stability`); a stability value needs at least two non-empty folds.
#' Runtime is the selection wall-clock, averaged over all folds.
#'
#' @param bm A `survfcbf_benchmark` from [nested_cv_select()].
#' @param theta Absolute-correlation similarity threshold, in `(0, 1]`. Default
#'   `0.9`.
#' @param seed Integer or `NULL`. Seed for the stability chance-correction sampling.
#' @param sim_mat Optional precomputed `p` by `p` absolute-correlation similarity
#'   matrix. `NULL` (default) computes it from the task according to `sim_scope`.
#'   Pass it in to avoid recomputing the matrix per selector (it is
#'   selector-independent); at `p` ~ 23k it is several GB, so recomputation is the
#'   dominant memory cost.
#' @param sim_scope How the similarity matrix is estimated when `sim_mat` is `NULL`.
#'   `"full"` (default) uses every row of the task; `"train"` uses the element-wise
#'   mean of the per-fold training-set correlation matrices, so the similarity
#'   structure never sees a held-out row. The matrix is unsupervised (it involves no
#'   outcome) and is used only to score agreement between already-selected sets, so
#'   `"full"` does not leak outcome information into a performance estimate, and a
#'   common reference across folds is what makes the per-fold sets comparable at
#'   all. `"train"` is the strict variant for studies where any use of held-out
#'   rows, supervised or not, must be ruled out; on data where the correlation
#'   structure is stable across folds the two agree. Ignored when `sim_mat` is
#'   given.
#'
#' @return A named list: `stability`, `theta`, `n_stability` (non-empty folds used),
#'   `parsimony_mean`, `parsimony_se`, `runtime_mean`, `runtime_se`.
#'
#' @seealso [nested_cv_select()], [score_benchmark()]
#' @export
benchmark_stability <- function(bm, theta = 0.9, seed = NULL, sim_mat = NULL,
                                sim_scope = c("full", "train")) {
  sim_scope <- match.arg(sim_scope)
  if (!requireNamespace("stabm", quietly = TRUE)) {
    stop("benchmark_stability requires the 'stabm' package.", call. = FALSE)
  }
  stopifnot("`bm` must be a survfcbf_benchmark" = inherits(bm, "survfcbf_benchmark"),
            "`theta` must be a single number in (0, 1]" =
              length(theta) == 1L && is.numeric(theta) && theta > 0 && theta <= 1)

  sets <- lapply(bm$folds, function(fd) as.integer(fd$selected))
  nonempty <- sets[vapply(sets, length, integer(1)) > 0L]
  stability <- if (length(nonempty) >= 2L) {
    if (is.null(sim_mat))                                   # selector-independent, so callers can precompute it
      sim_mat <- similarity_matrix(bm, scope = sim_scope)  # once and pass it in -- at p ~ 23k it is 4 GB and
    with_seed(seed, stabm::stabilityIntersectionCount(     # recomputing it per selector is what blows up memory.
      nonempty, sim.mat = sim_mat, threshold = theta))
  } else {
    NA_real_
  }

  card <- vapply(bm$folds, function(fd) as.numeric(fd$cardinality), numeric(1))
  card <- card[card > 0]                                   # parsimony over folds that selected something
  rt <- vapply(bm$folds, function(fd) as.numeric(fd$select_time), numeric(1))
  se <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_

  list(stability = as.numeric(stability), theta = theta, n_stability = length(nonempty),
       parsimony_mean = if (length(card)) mean(card) else NA_real_, parsimony_se = se(card),
       runtime_mean = mean(rt), runtime_se = se(rt))
}

#' Cluster-level recovery of a benchmark against the simulator truth
#'
#' For benchmarks on `simulate_dataset()` output, scores the per-fold selected sets
#' against the known cluster-level ground truth: mean sensitivity
#' (recovered driver clusters), false-discovery rate, and exact-recovery rate across
#' the outer folds. This is the "does it pick the right genes" evidence that real
#' data cannot give.
#'
#' @section Index spaces:
#' A fold's `selected` indices point into `task$feature_names`, which mlr3 returns
#' in *alphabetical* order (`V1`, `V10`, `V11`, ..., `V2`), whereas `manifest` rows
#' follow the original column order of `X`. Scoring `selected` against `manifest`
#' directly therefore reads the wrong rows. Recovery is scored by resolving each
#' fold's feature *names* against `feature_names`, so the two spaces cannot drift.
#'
#' @param bm A `survfcbf_benchmark` from [nested_cv_select()].
#' @param manifest The [make_sigma_blocks()] / `simulate_dataset()` manifest with
#'   `role` and `cluster_id`.
#' @param feature_names Character vector naming the columns of `X` in `manifest`
#'   row order. Defaults to `V1..Vp`, the names [surv_task()] generates for an
#'   unnamed simulator matrix. Supply the real `colnames(X)` when the task was
#'   built from a named matrix.
#'
#' @return A named list: `sensitivity_mean`, `sensitivity_se`, `fdr_mean`,
#'   `fdr_se`, `exact_recovery` (mean over folds), `n_recovery` (non-empty folds).
#'
#' @seealso [score_selection_clusterwise()], [benchmark_selectors()]
#' @export
benchmark_recovery <- function(bm, manifest, feature_names = NULL) {
  stopifnot("`bm` must be a survfcbf_benchmark" = inherits(bm, "survfcbf_benchmark"),
            "`manifest` must be a data.frame" = is.data.frame(manifest))
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(nrow(manifest)))
  stopifnot("`feature_names` must have one entry per `manifest` row" =
              is.character(feature_names) && length(feature_names) == nrow(manifest))
  sets <- lapply(bm$folds, function(fd) {
    if (is.null(fd$features)) {
      stop("fold has no `features`; recovery needs the selected names to resolve ",
           "the harness index space against `manifest` row order", call. = FALSE)
    }
    idx <- match(as.character(fd$features), feature_names)      # names -> manifest rows
    if (anyNA(idx)) {
      stop("selected features absent from `feature_names`: ",
           paste(utils::head(fd$features[is.na(idx)], 5L), collapse = ", "),
           call. = FALSE)
    }
    idx
  })
  sets <- sets[vapply(sets, length, integer(1)) > 0L]            # only folds that selected something
  se <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
  if (!length(sets)) {
    return(list(sensitivity_mean = NA_real_, sensitivity_se = NA_real_,
                fdr_mean = NA_real_, fdr_se = NA_real_, exact_recovery = NA_real_, n_recovery = 0L))
  }
  scs <- lapply(sets, function(s) score_selection_clusterwise(s, manifest))
  sens <- vapply(scs, function(x) x$sensitivity, numeric(1))
  fdr <- vapply(scs, function(x) x$fdr, numeric(1))
  exact <- vapply(scs, function(x) as.numeric(x$exact_recovery), numeric(1))
  list(sensitivity_mean = mean(sens), sensitivity_se = se(sens),
       fdr_mean = mean(fdr), fdr_se = se(fdr),
       exact_recovery = mean(exact), n_recovery = length(sets))
}

#' Standard selector set for the benchmark sweep
#'
#' Builds the named selector closures the benchmark compares:
#' survfcbf at its self-selected cardinality plus the nine [baseline_registry()]
#' baselines and, when a `learner` is supplied, the two inner-CV-tuned survfcbf
#' variants (`survfcbf_cvq`, with `q` tuned; `survfcbf_atprop`, with `q` and the
#' threshold tuned). When `topk` is given, a `survfcbf_topk` variant (the
#' explicit-cardinality mode for the low-event / high-dimensional regime) is added.
#' Ranking baselines are cut at `prop`; set-selecting ones keep their own cardinality.
#'
#' @param prop Cardinality proportion for the ranking baselines, in `(0, 1]`.
#'   Default `0.1`.
#' @param survfcbf_args Extra arguments to [survfcbf_select()] (e.g. `q`, `B`).
#' @param learner Optional `mlr3` survival learner. When given, the inner-CV-tuned
#'   survfcbf variants (via [survfcbf_cv_select()], tuned for this learner) are
#'   added; `NULL` (default) omits them.
#' @param q_grid,delta_grid Candidate grids for the tuned variants;
#'   `survfcbf_cvq` tunes `q_grid` only, `survfcbf_atprop` also tunes `delta_grid`.
#' @param cv_args Extra arguments to [survfcbf_cv_select()] (e.g. `inner_folds`, `B`).
#' @param topk Optional positive integer. When given, adds a `survfcbf_topk`
#'   selector: survfcbf in explicit-cardinality mode (`q = 2`, the validated
#'   genomic-regime discretisation) capped at `topk` predominant features. `NULL`
#'   (default) omits it.
#'
#' @return A named list of `(time, status, X) -> indices` selector closures.
#' @seealso [benchmark_selectors()], [baseline_registry()], [survfcbf_cv_select()]
#'
#' @examples
#' sels <- default_selectors(topk = 10)
#' names(sels)   # survfcbf in both cardinality modes, plus the baselines
#'
#' # each entry is a plain function, usable on its own
#' sim <- simulate_dataset(default_params(
#'   p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'   n_distractor_subblocks = 2L, distractor_size = 4L,
#'   calib_n_presample = 500, seed = 1))
#' sels$survfcbf_topk(sim$time, sim$status, sim$X)
#' @export
default_selectors <- function(prop = 0.1, survfcbf_args = list(), learner = NULL,
                              q_grid = c(4, 7, 10), delta_grid = c(0.02, 0.05, 0.1),
                              cv_args = list(), topk = NULL) {
  sels <- list(
    survfcbf = function(ti, st, X) do.call(survfcbf_select, c(list(ti, st, X), survfcbf_args))$selected
  )
  if (!is.null(topk)) {                             # explicit-cardinality survfcbf for the low-event regime
    sels[["survfcbf_topk"]] <- local({
      k <- as.integer(topk)
      function(ti, st, X) survfcbf_select(ti, st, X, q = 2L, top = k)$selected
    })
  }
  if (!is.null(learner)) {                          # learner-tuned survfcbf variants
    sels[["survfcbf_cvq"]] <- local({
      lr <- learner; qg <- q_grid; ca <- cv_args
      function(ti, st, X) do.call(survfcbf_cv_select, c(list(ti, st, X, lr, q_grid = qg), ca))$selected
    })
    sels[["survfcbf_atprop"]] <- local({
      lr <- learner; qg <- q_grid; dg <- delta_grid; ca <- cv_args
      function(ti, st, X) do.call(survfcbf_cv_select, c(list(ti, st, X, lr, q_grid = qg, delta_grid = dg), ca))$selected
    })
  }
  for (e in baseline_registry()) {
    sels[[e$name]] <- local({
      entry <- e
      pa <- if (e$self_selects) NULL else prop
      function(ti, st, X) run_baseline(entry, ti, st, X, prop = pa)$selected
    })
  }
  sels
}

#' Sweep selectors through the nested-CV benchmark
#'
#' Runs each selector through [nested_cv_select()], scores it with
#' [score_benchmark()] and [benchmark_stability()], and assembles one row per
#' selector: the discrimination and calibration metrics (with MC
#' standard errors), stability, parsimony, runtime and the failed-fold count.
#'
#' @param task An `mlr3proba::TaskSurv` (see [surv_task()]).
#' @param selectors A named list of `(time, status, X) -> indices` closures, e.g.
#'   from [default_selectors()].
#' @param learner An `mlr3` survival learner.
#' @param outer_folds Outer CV folds. Default `10`.
#' @param measures Measures for [score_benchmark()]. `NULL` uses the default set.
#' @param truncation Passed to [score_benchmark()]: the training-time quantile at
#'   which the truncated Uno's C (`uno_c_trunc`, the primary discrimination metric)
#'   is cut. Default `0.8`, `NULL` to skip it.
#' @param theta Stability similarity threshold. Default `0.9`.
#' @param manifest Optional `simulate_dataset()` manifest. When given, cluster-level
#'   recovery columns (`sensitivity`, `fdr`, `exact_recovery`) are added from
#'   [benchmark_recovery()] (the known-truth evidence on simulated data).
#' @param seed Integer or `NULL`. Master seed, shared across selectors so they see
#'   the same folds.
#' @param sim_scope How the stability similarity matrix is estimated, `"full"`
#'   (default) or `"train"`; see [similarity_matrix()]. The matrix is built once
#'   and shared across selectors, except for `"train"` with a `NULL` seed, where
#'   each selector draws its own folds and the matrix is rebuilt per selector.
#' @param rng Per-fold seed derivation, `"sequential"` (default) or `"stream"`;
#'   passed to [nested_cv_select()].
#'
#' @return A data.frame with one row per selector: `selector`, the per-measure mean
#'   and `_se` columns, `stability`, `parsimony` (`+_se`), `runtime` (`+_se`),
#'   `n_failed`, `n_null_fit` (folds whose learner fell back to the null model, see
#'   [score_benchmark()]), and (with `manifest`) `sensitivity`, `fdr`,
#'   `exact_recovery`.
#'
#' @param feature_names Passed to [benchmark_recovery()]: the columns of `X` in
#'   `manifest` row order. Only consulted when `manifest` is given.
#'
#' @seealso [default_selectors()], [nested_cv_select()], [benchmark_recovery()]
#'
#' @examples
#' if (requireNamespace("mlr3proba", quietly = TRUE) &&
#'     requireNamespace("stabm", quietly = TRUE)) {
#'   sim <- simulate_dataset(default_params(
#'     p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'     n_distractor_subblocks = 2L, distractor_size = 4L,
#'     calib_n_presample = 500, seed = 1))
#'   task <- surv_task(sim$time, sim$status, sim$X)
#'   sels <- default_selectors(topk = 5)[c("survfcbf_topk", "variance")]
#'
#'   # with the simulator's manifest, cluster-level recovery is scored too
#'   res <- benchmark_selectors(task, sels, mlr3::lrn("surv.coxph"),
#'                              outer_folds = 3, manifest = sim$manifest, seed = 1)
#'   res[, c("selector", "uno_c_trunc", "stability", "parsimony", "sensitivity", "fdr")]
#' }
#' @export
benchmark_selectors <- function(task, selectors, learner, outer_folds = 10L,
                                measures = NULL, truncation = 0.8, theta = 0.9,
                                manifest = NULL, seed = NULL,
                                sim_scope = c("full", "train"), rng = c("sequential", "stream"),
                                feature_names = NULL) {
  sim_scope <- match.arg(sim_scope)
  rng <- match.arg(rng)
  stopifnot("`selectors` must be a named list of functions" =
              is.list(selectors) && length(selectors) >= 1L && !is.null(names(selectors)))

  ## Feature-similarity matrix for the stability metric, built lazily from the first
  ## benchmark and shared across selectors (it is selector-independent, and at p ~ 23k
  ## recomputing it per selector is a 4 GB-per-call OOM). Sharing is only valid for
  ## "train" scope when a seed fixes the folds, since the matrix then depends on the
  ## fold assignment; with a NULL seed each selector draws its own folds, so the
  ## train-scoped matrix is rebuilt per selector rather than shared incorrectly.
  sim_cache <- new.env(parent = emptyenv())
  sim_for <- function(bm) {
    if (sim_scope == "train" && is.null(seed)) return(similarity_matrix(bm, scope = "train"))
    if (is.null(sim_cache$mat)) sim_cache$mat <- similarity_matrix(bm, scope = sim_scope)
    sim_cache$mat
  }

  rows <- lapply(names(selectors), function(nm) {
    bm <- nested_cv_select(task, selectors[[nm]], learner, outer_folds = outer_folds,
                           seed = seed, rng = rng)
    sc <- score_benchmark(bm, measures = measures, truncation = truncation)
    st <- benchmark_stability(bm, theta = theta, seed = seed, sim_mat = sim_for(bm))
    metric_cols <- do.call(data.frame, c(
      as.list(stats::setNames(sc$summary$mean, sc$summary$metric)),
      stats::setNames(as.list(sc$summary$se), paste0(sc$summary$metric, "_se")),
      list(check.names = FALSE)
    ))
    row <- data.frame(selector = nm, metric_cols,
                      stability = st$stability,
                      parsimony = st$parsimony_mean, parsimony_se = st$parsimony_se,
                      runtime = st$runtime_mean, runtime_se = st$runtime_se,
                      n_failed = sc$n_failed, n_null_fit = sc$n_null_fit,
                      row.names = NULL, check.names = FALSE)
    if (!is.null(manifest)) {
      rc <- benchmark_recovery(bm, manifest, feature_names = feature_names)
      row$sensitivity <- rc$sensitivity_mean
      row$fdr <- rc$fdr_mean
      row$exact_recovery <- rc$exact_recovery
    }
    row
  })
  do.call(rbind, rows)
}
