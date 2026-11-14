# Acceptance tests for the two-stage permutation-assisted group-Lasso ensemble.
# Stage 1 (WRF ensemble of base selectors over folds) + stage 2 (PAGL) is costly, so one
# shared simulation and small cv_folds / K keep it lean.

gl_sim <- simulate_dataset(default_params(
  p = 120, k_drivers = 3, cluster_size = 4, n = 500,
  n_distractor_subblocks = 4, distractor_size = 6,
  censor_fraction = 0.10, target_C = 0.90, seed = 31, calib_n_presample = 1000))

test_that("the two-stage ensemble recovers signal and groups correlated drivers", {
  skip_if_not_installed("grpreg")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("praznik")
  r <- grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, cv_folds = 3, K = 8, seed = 1)

  expect_length(r$score, 120L)
  expect_length(r$group, 120L)
  expect_true(length(r$candidate) >= length(r$selected))   # stage-1 candidates contain the final set
  expect_true(all(r$selected %in% r$candidate))
  # the correlated driver clusters collapse into a few groups among the candidates
  drv <- which(gl_sim$manifest$role %in% c("driver", "redundant"))
  expect_lte(length(setdiff(unique(r$group[drv]), 0L)), 5L)
  sc <- score_selection_clusterwise(r$selected, gl_sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
  expect_lte(sc$fdr, 0.3)
})

test_that("the ensemble cuts at prop and is reproducible under a seed", {
  skip_if_not_installed("grpreg")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("praznik")
  a <- grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, prop = 0.1, cv_folds = 2, K = 4, seed = 2)
  b <- grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, prop = 0.1, cv_folds = 2, K = 4, seed = 2)
  expect_identical(a$selected, b$selected)
  expect_length(a$selected, max(1L, round(0.1 * 120)))
})

test_that("grplasso_ensemble validates inputs before fitting", {
  skip_if_not_installed("grpreg")
  expect_error(grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, rho_t = 1.5), "rho_t")
  expect_error(grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, cv_folds = 1), "cv_folds")
  expect_error(grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, ensemble = 3), "ensemble")
  expect_error(grplasso_ensemble_select(gl_sim$time, gl_sim$status, gl_sim$X, tau = 0), "tau")
})
