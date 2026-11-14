# Acceptance tests for the embedded selectors (cox-glmnet, RSF).
# One shared simulation (built once) keeps the RSF permutation-importance cost down.

emb_sim <- simulate_dataset(default_params(
  p = 150, k_drivers = 3, cluster_size = 4, n = 500,
  n_distractor_subblocks = 5, distractor_size = 8,
  censor_fraction = 0.10, target_C = 0.90, seed = 31, calib_n_presample = 1000))
emb_sc <- function(s) score_selection_clusterwise(s, emb_sim$manifest)

test_that("cox_glmnet lasso self-selects a nonzero set and recovers signal", {
  skip_if_not_installed("glmnet")
  r <- cox_glmnet_select(emb_sim$time, emb_sim$status, emb_sim$X, alpha = 1, seed = 1)

  expect_length(r$score, 150L)
  expect_equal(r$score, abs(r$coef))
  expect_setequal(r$selected, which(r$coef != 0))     # lasso: the nonzero coefficients are the set
  expect_gte(length(r$selected), 3L)
  expect_gte(emb_sc(r$selected)$sensitivity, 2 / 3)   # recovers clusters (FDR is the known lasso weakness)
})

test_that("cox_glmnet ridge is a pure ranker, cut at prop", {
  skip_if_not_installed("glmnet")
  ridge <- cox_glmnet_select(emb_sim$time, emb_sim$status, emb_sim$X, alpha = 0, seed = 1)
  expect_null(ridge$selected)                         # ridge has no intrinsic sparsity
  rp <- cox_glmnet_select(emb_sim$time, emb_sim$status, emb_sim$X, alpha = 0, prop = 0.05, seed = 1)
  expect_length(rp$selected, max(1L, round(0.05 * 150)))
  expect_gte(emb_sc(rp$selected)$sensitivity, 2 / 3)
})

test_that("rsf_select ranks by permutation importance and recovers signal", {
  skip_if_not_installed("ranger")
  r <- rsf_select(emb_sim$time, emb_sim$status, emb_sim$X, prop = 0.05, num.trees = 100, seed = 1)
  expect_length(r$score, 150L)
  expect_null(rsf_select(emb_sim$time, emb_sim$status, emb_sim$X, num.trees = 50, seed = 1)$selected)
  sc <- emb_sc(r$selected)
  expect_gte(sc$sensitivity, 2 / 3)
  expect_lte(sc$fdr, 0.4)
})

test_that("embedded selectors are reproducible under a fixed seed", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ranger")
  expect_identical(
    cox_glmnet_select(emb_sim$time, emb_sim$status, emb_sim$X, alpha = 1, seed = 7)$selected,
    cox_glmnet_select(emb_sim$time, emb_sim$status, emb_sim$X, alpha = 1, seed = 7)$selected)
  expect_identical(
    rsf_select(emb_sim$time, emb_sim$status, emb_sim$X, prop = 0.05, num.trees = 50, seed = 7)$score,
    rsf_select(emb_sim$time, emb_sim$status, emb_sim$X, prop = 0.05, num.trees = 50, seed = 7)$score)
})

test_that("cox_glmnet validates alpha and shared dimensions", {
  skip_if_not_installed("glmnet")
  expect_error(cox_glmnet_select(emb_sim$time, emb_sim$status, emb_sim$X, alpha = 2), "alpha")
  expect_error(cox_glmnet_select(emb_sim$time[-1], emb_sim$status, emb_sim$X), "share nrow")
})
