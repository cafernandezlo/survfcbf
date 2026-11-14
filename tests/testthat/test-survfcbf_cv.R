# Acceptance tests for the inner-CV-tuned survfcbf selector.

if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")

cv_sim <- simulate_dataset(default_params(
  p = 60, n = 400, k_drivers = 3, cluster_size = 4,
  n_distractor_subblocks = 3, distractor_size = 4, target_C = 0.90,
  censor_fraction = 0.20, seed = 7, calib_n_presample = 800))

test_that("survfcbf_cv_select tunes q by inner CV and recovers signal", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  r <- survfcbf_cv_select(cv_sim$time, cv_sim$status, cv_sim$X, lrn, q_grid = c(5, 10),
                          inner_folds = 3, B = 30, seed = 1)

  expect_true(r$q %in% c(5, 10))                       # chose a q from the grid
  expect_gte(length(r$selected), 1L)
  expect_equal(nrow(r$cv), 2L)
  expect_true(all(c("q", "delta", "score") %in% names(r$cv)))
  sc <- score_selection_clusterwise(r$selected, cv_sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
})

test_that("a delta_grid tunes the threshold (the at-prop variant)", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  r <- survfcbf_cv_select(cv_sim$time, cv_sim$status, cv_sim$X, lrn, q_grid = 10,
                          delta_grid = c(0.02, 0.05, 0.1), inner_folds = 3, B = 30, seed = 1)

  expect_true(r$delta %in% c(0.02, 0.05, 0.1))         # chose a delta from the grid
  expect_gte(length(r$selected), 1L)
  expect_equal(nrow(r$cv), 3L)
})

test_that("survfcbf_cv_select is reproducible and validates inputs", {
  skip_if_not_installed("mlr3proba")
  skip_if_not_installed("mlr3extralearners")
  skip_if_not_installed("glmnet")
  lrn <- mlr3::lrn("surv.cv_glmnet", alpha = 0)
  a <- survfcbf_cv_select(cv_sim$time, cv_sim$status, cv_sim$X, lrn, q_grid = 5, inner_folds = 2, B = 20, seed = 2)
  b <- survfcbf_cv_select(cv_sim$time, cv_sim$status, cv_sim$X, lrn, q_grid = 5, inner_folds = 2, B = 20, seed = 2)
  expect_identical(a$selected, b$selected)
  expect_equal(a$q, b$q)
  expect_error(survfcbf_cv_select(cv_sim$time, cv_sim$status, cv_sim$X, lrn, inner_folds = 1), "inner_folds")
})
