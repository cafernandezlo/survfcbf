# Acceptance tests for the leakage-free nested-CV harness.

if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")

test_that("surv_task builds a TaskSurv from the matrix form", {
  skip_if_not_installed("mlr3proba")
  sim <- simulate_dataset(default_params(
    p = 40, n = 200, seed = 1, calib_n_presample = 600,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4))
  task <- surv_task(sim$time, sim$status, sim$X)
  expect_s3_class(task, "TaskSurv")
  expect_equal(task$nrow, 200L)
  expect_equal(length(task$feature_names), 40L)
  expect_setequal(task$target_names, c("time", "status"))
})

test_that("nested_cv_select runs leakage-free and carries signal end to end", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  sim <- simulate_dataset(default_params(
    p = 60, n = 300, seed = 2, calib_n_presample = 800,
    k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 3, distractor_size = 4,
    target_C = 0.90, censor_fraction = 0.10))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)                       # L2-Cox downstream
  sfn <- function(ti, st, X) cox_glmnet_select(ti, st, X, alpha = 1)$selected   # lasso, supervised, stable

  bm <- nested_cv_select(task, sfn, lrn, outer_folds = 3, seed = 1)
  expect_s3_class(bm, "survfcbf_benchmark")
  expect_length(bm$folds, 3L)

  # fold mechanics hold for every fold (selection and partition do not depend on the learner)
  for (fd in bm$folds) {
    expect_true(all(fd$selected >= 1L & fd$selected <= 60L))
    expect_length(intersect(fd$train_set, fd$test_set), 0L)          # train/test disjoint -> no leakage
    expect_gte(fd$cardinality, 1L)
  }
  expect_setequal(sort(unlist(lapply(bm$folds, `[[`, "test_set"))), task$row_ids)  # folds partition rows

  # successful folds carry signal (mean Harrell C > 0.5)
  ok <- Filter(function(fd) !is.null(fd$prediction), bm$folds)
  expect_gte(length(ok), 2L)
  hc <- mean(vapply(ok, function(fd) fd$prediction$score(mlr3::msr("surv.cindex")), numeric(1)))
  expect_gt(hc, 0.5)
})

test_that("an empty selection is recorded as a failed fold, not masked as no-filter", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  sim <- simulate_dataset(default_params(
    p = 30, n = 150, seed = 4, calib_n_presample = 500,
    k_drivers = 2, cluster_size = 2, n_distractor_subblocks = 2, distractor_size = 3))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)

  bm <- nested_cv_select(task, function(ti, st, X) integer(0), lrn, outer_folds = 3, seed = 1)
  for (fd in bm$folds) {
    expect_equal(fd$cardinality, 0L)              # the selector chose nothing; not inflated to p
    expect_null(fd$prediction)
    expect_length(fd$selected, 0L)
    expect_match(fd$error, "empty")
  }
})

test_that("default_surv_measures provides Uno, Harrell and IBS", {
  skip_if_not_installed("mlr3proba")
  ms <- default_surv_measures()
  expect_setequal(names(ms), c("uno_c", "harrell_c", "ibs"))
  expect_true(all(vapply(ms, function(m) inherits(m, "Measure"), logical(1))))
})

test_that("score_benchmark aggregates fold metrics with MC standard errors", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  sim <- simulate_dataset(default_params(
    p = 50, n = 300, seed = 2, calib_n_presample = 700,
    k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 3, distractor_size = 3,
    target_C = 0.90, censor_fraction = 0.10))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  bm <- nested_cv_select(task, function(ti, st, X) cox_glmnet_select(ti, st, X, alpha = 1)$selected,
                         lrn, outer_folds = 3, seed = 1)
  sc <- score_benchmark(bm)

  expect_setequal(sc$summary$metric, c("uno_c", "harrell_c", "ibs", "uno_c_trunc"))
  expect_equal(nrow(sc$per_fold), 3L)
  expect_true(all(c("uno_c", "harrell_c", "ibs", "uno_c_trunc", "t_max", "cardinality",
                    "select_time") %in% names(sc$per_fold)))
  expect_gt(sc$summary$mean[sc$summary$metric == "uno_c"], 0.5)       # signal recovered
  expect_gt(sc$summary$mean[sc$summary$metric == "uno_c_trunc"], 0.5)
  expect_lt(sc$summary$mean[sc$summary$metric == "ibs"], 0.3)
  expect_true(all(is.finite(sc$summary$se)))
  expect_equal(sc$n_failed, 0L)
})

test_that("score_benchmark truncates Uno's C at a horizon learned on each training fold", {
  skip_if_not_installed("mlr3proba")
  sim <- simulate_dataset(default_params(
    p = 30, n = 240, seed = 3, calib_n_presample = 600,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 3,
    censor_fraction = 0.5))
  task <- surv_task(sim$time, sim$status, sim$X)
  bm <- nested_cv_select(task, function(ti, st, X) survfcbf_select(ti, st, X, q = 2, top = 3)$selected,
                         mlr3::lrn("surv.coxph"), outer_folds = 3, seed = 1)
  sc <- score_benchmark(bm, truncation = 0.7)

  ## the horizon is the quantile of the training rows' observed times, never the test rows'
  for (f in seq_along(bm$folds)) {
    fd <- bm$folds[[f]]
    expect_equal(sc$per_fold$t_max[f], unname(stats::quantile(sim$time[fd$train_set], 0.7)))
    m <- mlr3::msr("surv.cindex", weight_meth = "G2", t_max = sc$per_fold$t_max[f])
    expect_equal(sc$per_fold$uno_c_trunc[f],
                 unname(fd$prediction$score(m, task = bm$task, train_set = fd$train_set)))
  }
  expect_equal(sc$summary$mean[sc$summary$metric == "uno_c_trunc"], mean(sc$per_fold$uno_c_trunc))

  ## truncation = NULL gives the untruncated set alone, with the same values
  sc0 <- score_benchmark(bm, truncation = NULL)
  expect_setequal(sc0$summary$metric, c("uno_c", "harrell_c", "ibs"))
  expect_false(any(c("uno_c_trunc", "t_max") %in% names(sc0$per_fold)))
  expect_equal(sc0$per_fold$uno_c, sc$per_fold$uno_c)

  expect_error(score_benchmark(bm, truncation = 0), "truncation")
  expect_error(score_benchmark(bm, truncation = c(0.5, 0.8)), "truncation")
  clash <- list(uno_c_trunc = mlr3::msr("surv.cindex", id = "uno_c_trunc"))
  expect_error(score_benchmark(bm, measures = clash), "uno_c_trunc")
})

test_that("score_benchmark records failed folds as NA, not hidden", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  sim <- simulate_dataset(default_params(
    p = 30, n = 150, seed = 4, calib_n_presample = 500,
    k_drivers = 2, cluster_size = 2, n_distractor_subblocks = 2, distractor_size = 3))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  bm <- nested_cv_select(task, function(ti, st, X) integer(0), lrn, outer_folds = 3, seed = 1)
  sc <- score_benchmark(bm)
  expect_equal(sc$n_failed, 3L)
  expect_true(all(sc$summary$n_folds == 0L))
  expect_true(all(is.na(sc$summary$mean)))
  expect_true(all(is.na(sc$per_fold$null_fit)))   # no prediction -> no null-fit verdict
  expect_equal(sc$n_null_fit, 0L)
})

# Low-event fixture for the null-fit tests (~60 events in 200, p = 200): the regime
# where the 1se rule sends the ridge Cox to the first lambda of its path. Two panels
# of the same task: the true drivers and twenty pure-noise features.
null_fit_fixture <- function() {
  sim <- simulate_dataset(default_params(
    p = 200, n = 200, seed = 1, calib_n_presample = 600,
    k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 3, distractor_size = 4,
    censor_fraction = 0.7))
  task <- surv_task(sim$time, sim$status, sim$X)
  ix <- function(role) match(paste0("V", which(sim$manifest$role %in% role)), task$feature_names)
  list(task = task, drivers = ix(c("driver", "interaction")), noise = ix("noise")[1:20])
}

test_that("score_benchmark flags folds whose learner fell back to the null model", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  fx <- null_fit_fixture()
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)                  # default s = "lambda.1se"
  sc_d <- score_benchmark(nested_cv_select(fx$task, function(ti, st, X) fx$drivers, lrn,
                                           outer_folds = 3, seed = 1))
  sc_n <- score_benchmark(nested_cv_select(fx$task, function(ti, st, X) fx$noise, lrn,
                                           outer_folds = 3, seed = 1))
  expect_true(all(sc_d$per_fold$null_fit))                       # every fold is the null model
  expect_equal(sc_d$n_null_fit, 3L)
  expect_equal(sc_n$n_null_fit, 3L)
  expect_identical(sc_d$per_fold$ibs, sc_n$per_fold$ibs)         # the symptom: IBS ignores the panel
  expect_false(identical(sc_d$per_fold$uno_c, sc_n$per_fold$uno_c))   # while the C-index still moves
})

test_that("with lambda.min the ridge Cox distribution depends on the selected panel", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  fx <- null_fit_fixture()
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")
  bm_d <- nested_cv_select(fx$task, function(ti, st, X) fx$drivers, lrn, outer_folds = 3, seed = 1)
  sc_d <- score_benchmark(bm_d)
  sc_n <- score_benchmark(nested_cv_select(fx$task, function(ti, st, X) fx$noise, lrn,
                                           outer_folds = 3, seed = 1))
  expect_equal(sc_d$n_null_fit, 0L)                              # the drivers yield a real model
  for (fd in bm_d$folds) {                                       # curves vary across held-out rows
    S <- fd$prediction$data$distr
    expect_gt(max(abs(S - rep(S[1L, ], each = nrow(S)))), 1e-3)
  }
  expect_false(isTRUE(all.equal(sc_d$per_fold$ibs, sc_n$per_fold$ibs)))
  expect_lt(mean(sc_d$per_fold$ibs), mean(sc_n$per_fold$ibs))   # the true drivers calibrate better
})

test_that("benchmark_stability computes SMA-Count stability, parsimony and runtime", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("stabm")
  sim <- simulate_dataset(default_params(
    p = 40, n = 200, seed = 5, calib_n_presample = 500,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4))
  task <- surv_task(sim$time, sim$status, sim$X)
  # mock a benchmark of three overlapping selections to exercise the stability layer fast
  mk <- function(sel, t) list(selected = sel, cardinality = length(sel), select_time = t)
  bm <- structure(list(task = task, outer_folds = 3, fold_id = NULL,
                       folds = list(mk(c(1, 2, 3, 4), 0.10), mk(c(1, 2, 3, 5), 0.12), mk(c(1, 2, 4, 6), 0.08))),
                  class = "survfcbf_benchmark")

  st <- benchmark_stability(bm, theta = 0.9, seed = 1)
  expect_true(is.finite(st$stability) && st$stability >= -1 && st$stability <= 1)
  expect_equal(st$n_stability, 3L)
  expect_equal(st$parsimony_mean, 4)
  expect_gt(st$runtime_mean, 0)
  expect_equal(st$stability, benchmark_stability(bm, theta = 0.9, seed = 1)$stability)   # reproducible
})

test_that("rng = 'stream' gives independent, reproducible per-fold streams", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  sim <- simulate_dataset(default_params(
    p = 40, n = 200, seed = 5, calib_n_presample = 500,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  stoch <- function(ti, st, X) sort(sample.int(ncol(X), 5L))   # stochastic: exercises the RNG path

  # reproducible: identical runs give identical selections
  b1 <- nested_cv_select(task, stoch, lrn, outer_folds = 3, seed = 7, rng = "stream")
  b2 <- nested_cv_select(task, stoch, lrn, outer_folds = 3, seed = 7, rng = "stream")
  expect_identical(lapply(b1$folds, `[[`, "selected"), lapply(b2$folds, `[[`, "selected"))

  # the fold assignment is the same under both modes (only the draws differ)
  b3 <- nested_cv_select(task, stoch, lrn, outer_folds = 3, seed = 7, rng = "sequential")
  expect_identical(b1$fold_id, b3$fold_id)

  # a stream needs a master seed
  expect_error(nested_cv_select(task, stoch, lrn, outer_folds = 3, rng = "stream"), "seed")

  # the global RNG kind AND state are restored after a stream run, so a stream benchmark
  # does not perturb the caller's random sequence (restoring the kind re-inits the
  # generator, so the seed has to be put back *after* the kind, not before)
  RNGkind("Mersenne-Twister"); set.seed(1)
  before <- .Random.seed
  expected <- runif(3)
  RNGkind("Mersenne-Twister"); set.seed(1)
  invisible(nested_cv_select(task, stoch, lrn, outer_folds = 3, seed = 7, rng = "stream"))
  expect_equal(RNGkind()[1], "Mersenne-Twister")
  expect_identical(.Random.seed, before)
  expect_equal(runif(3), expected)
})

test_that("similarity_matrix supports full and train scopes and stability accepts both", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("stabm")
  sim <- simulate_dataset(default_params(
    p = 40, n = 200, seed = 5, calib_n_presample = 500,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4))
  task <- surv_task(sim$time, sim$status, sim$X)
  mk <- function(sel, t) list(selected = sel, cardinality = length(sel), select_time = t)
  bm <- structure(list(task = task, outer_folds = 3, fold_id = rep(1:3, length.out = 200),
                       folds = list(mk(c(1, 2, 3, 4), 0.1), mk(c(1, 2, 3, 5), 0.1), mk(c(1, 2, 4, 6), 0.1))),
                  class = "survfcbf_benchmark")

  Sfull  <- similarity_matrix(bm, "full")
  Strain <- similarity_matrix(bm, "train")
  expect_equal(dim(Sfull), c(40L, 40L))
  expect_equal(dim(Strain), c(40L, 40L))
  expect_equal(unname(diag(Strain)), rep(1, 40))
  expect_true(all(Strain >= 0 & Strain <= 1))
  expect_true(isSymmetric(unname(Strain)))
  # each training fold keeps ~2/3 of the rows, so the train-scoped estimate stays close
  expect_lt(mean(abs(Strain - Sfull)), 0.1)
  # deterministic: same fold assignment, same matrix
  expect_identical(Strain, similarity_matrix(bm, "train"))

  st_train <- benchmark_stability(bm, theta = 0.9, seed = 1, sim_scope = "train")
  expect_true(is.finite(st_train$stability) && st_train$stability >= -1 && st_train$stability <= 1)

  # a benchmark without a fold assignment cannot be train-scoped, and says so
  bm$fold_id <- NULL
  expect_error(similarity_matrix(bm, "train"), "fold assignment")
})

test_that("benchmark_selectors passes sim_scope through", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("stabm")
  skip_if_not_installed("glmnet")
  sim <- simulate_dataset(default_params(
    p = 40, n = 200, seed = 5, calib_n_presample = 500,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  sels <- list(variance = function(ti, st, X) cut_ranking(apply(X, 2L, stats::var), k = 5L))
  tab <- benchmark_selectors(task, sels, lrn, outer_folds = 3L, seed = 1L, sim_scope = "train")
  expect_equal(nrow(tab), 1L)
  expect_true(is.finite(tab$stability))
})

test_that("benchmark_stability returns NA stability when folds select nothing", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("stabm")
  sim <- simulate_dataset(default_params(
    p = 30, n = 150, seed = 6, calib_n_presample = 400,
    k_drivers = 2, cluster_size = 2, n_distractor_subblocks = 2, distractor_size = 3))
  task <- surv_task(sim$time, sim$status, sim$X)
  mk0 <- function() list(selected = integer(0), cardinality = 0L, select_time = 0.05)
  bm <- structure(list(task = task, outer_folds = 3, fold_id = NULL,
                       folds = list(mk0(), mk0(), mk0())), class = "survfcbf_benchmark")
  st <- benchmark_stability(bm, seed = 1)
  expect_true(is.na(st$stability))
  expect_equal(st$n_stability, 0L)
  expect_true(is.na(st$parsimony_mean))
})

test_that("benchmark_recovery scores cluster-level sensitivity and FDR vs the manifest", {
  skip_if_not_installed("mlr3proba")
  sim <- simulate_dataset(default_params(
    p = 60, n = 300, seed = 7, calib_n_presample = 700,
    k_drivers = 3, cluster_size = 4, n_distractor_subblocks = 3, distractor_size = 4,
    target_C = 0.90))
  task <- surv_task(sim$time, sim$status, sim$X)
  drv <- which(sim$manifest$role == "driver")             # one member per driver cluster
  nm <- paste0("V", seq_len(nrow(sim$manifest)))          # manifest row order
  mk <- function(sel) list(selected = sel, features = nm[sel],
                           cardinality = length(sel), select_time = 0.1)
  bm <- structure(list(task = task, outer_folds = 3, fold_id = NULL,
                       folds = list(mk(drv), mk(drv), mk(drv))), class = "survfcbf_benchmark")

  rc <- benchmark_recovery(bm, sim$manifest)
  expect_equal(rc$sensitivity_mean, 1)                    # all driver clusters recovered
  expect_equal(rc$fdr_mean, 0)                            # only drivers selected
  expect_equal(rc$n_recovery, 3L)
})

test_that("benchmark_recovery resolves the harness index space, not raw positions", {
  skip_if_not_installed("mlr3proba")
  # Regression: mlr3 returns feature_names alphabetically (V1, V10, V11, ..., V2), so a
  # fold's `selected` indices are NOT manifest row numbers. Scoring them directly made an
  # oracle selector look like it had picked noise. Route through nested_cv_select so the
  # reordering is actually exercised.
  sim <- simulate_dataset(default_params(
    p = 30, n = 200, seed = 3, calib_n_presample = 500,
    k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4))
  task <- surv_task(sim$time, sim$status, sim$X)
  nm <- paste0("V", seq_len(nrow(sim$manifest)))
  drv_names <- nm[sim$manifest$role == "driver"]
  expect_false(identical(task$feature_names, nm))         # the reordering is real

  oracle <- function(ti, st, X) which(colnames(X) %in% drv_names)
  bm <- nested_cv_select(task, oracle, mlr3::lrn("surv.kaplan"), outer_folds = 3, seed = 1)
  expect_setequal(bm$folds[[1]]$features, drv_names)      # the oracle did pick the drivers

  rc <- benchmark_recovery(bm, sim$manifest)
  expect_equal(rc$sensitivity_mean, 1)                    # ... and recovery must say so
  expect_equal(rc$fdr_mean, 0)
  expect_equal(rc$exact_recovery, 1)
})

test_that("default_selectors builds survfcbf plus the registry baselines", {
  skip_if_not_installed("mlr3proba")
  sels <- default_selectors()
  expect_length(sels, 10L)
  expect_true(all(c("survfcbf", "variance", "carss", "rsf", "grplasso_ensemble") %in% names(sels)))
  expect_true(all(vapply(sels, is.function, logical(1))))
})

test_that("default_selectors adds the tuned survfcbf variants when given a learner", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  sels <- default_selectors(learner = mlr3::lrn("surv.cv_glmnet", alpha = 0))
  expect_length(sels, 12L)
  expect_true(all(c("survfcbf", "survfcbf_cvq", "survfcbf_atprop") %in% names(sels)))
  expect_true(all(vapply(sels, is.function, logical(1))))
})

test_that("default_selectors adds a survfcbf_topk variant capped at topk", {
  sels <- default_selectors(topk = 5)
  expect_length(sels, 11L)                                  # survfcbf + survfcbf_topk + 9 baselines
  expect_true("survfcbf_topk" %in% names(sels))
  # the closure runs in explicit-cardinality mode and respects the cap
  sim <- simulate_dataset(default_params(
    p = 150, k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 2,
    distractor_size = 4, n = 250, censor_fraction = 0.6, seed = 9, calib_n_presample = 800))
  idx <- sels[["survfcbf_topk"]](sim$time, sim$status, sim$X)
  expect_lte(length(idx), 5L)
  expect_gte(length(idx), 1L)
})

test_that("benchmark_selectors assembles a one-row-per-selector table", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("stabm")
  sim <- simulate_dataset(default_params(
    p = 50, n = 300, seed = 2, calib_n_presample = 700,
    k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 3, distractor_size = 3,
    target_C = 0.90, censor_fraction = 0.10))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  sels <- list(
    variance = function(ti, st, X) variance_select(ti, st, X, prop = 0.1)$selected,
    cox_l1   = function(ti, st, X) cox_glmnet_select(ti, st, X, alpha = 1)$selected
  )
  tab <- benchmark_selectors(task, sels, lrn, outer_folds = 3, manifest = sim$manifest, seed = 1)

  expect_equal(nrow(tab), 2L)
  expect_true(all(c("selector", "uno_c", "harrell_c", "ibs", "stability",
                    "parsimony", "runtime", "n_failed", "n_null_fit") %in% names(tab)))
  expect_true(all(c("sensitivity", "fdr", "exact_recovery") %in% names(tab)))   # manifest -> recovery cols
  expect_true(all(tab$n_null_fit >= 0L & tab$n_null_fit <= 3L))
  expect_setequal(tab$selector, c("variance", "cox_l1"))
  # the variance filter is unsupervised (Uno ~ 0.5), the supervised lasso discriminates
  expect_lt(abs(tab$uno_c[tab$selector == "variance"] - 0.5), 0.1)
  expect_gt(tab$uno_c[tab$selector == "cox_l1"], 0.6)
  expect_true(all(is.finite(tab$stability)))
  expect_true(all(is.finite(tab$sensitivity)))
})

test_that("nested_cv_select is reproducible under a fixed seed", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  sim <- simulate_dataset(default_params(
    p = 40, n = 200, seed = 3, calib_n_presample = 600,
    k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 4,
    target_C = 0.90))
  task <- surv_task(sim$time, sim$status, sim$X)
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  sfn <- function(ti, st, X) survfcbf_select(ti, st, X, q = 8, B = 10)$selected   # stochastic selection

  a <- nested_cv_select(task, sfn, lrn, outer_folds = 3, seed = 7)
  b <- nested_cv_select(task, sfn, lrn, outer_folds = 3, seed = 7)
  expect_identical(a$fold_id, b$fold_id)
  expect_identical(lapply(a$folds, `[[`, "selected"), lapply(b$folds, `[[`, "selected"))
})
