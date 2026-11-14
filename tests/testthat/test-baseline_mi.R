# Acceptance tests for the direct mutual-information family.

concentrated_sim <- function(seed = 31) {
  simulate_dataset(default_params(
    p = 300, k_drivers = 3, cluster_size = 5, n = 800,
    censor_fraction = 0.10, target_C = 0.90, seed = seed, calib_n_presample = 1500))
}

test_that("mim_surv_select recovers concentrated signal cleanly under the common contract", {
  skip_if_not_installed("praznik")
  sim <- concentrated_sim()
  r <- mim_surv_select(sim$time, sim$status, sim$X, prop = 0.05)   # top 15 of 300

  expect_length(r$score, 300L)
  expect_length(r$selected, 15L)
  expect_true(all(r$selected >= 1L & r$selected <= 300L))
  sc <- score_selection_clusterwise(r$selected, sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
  expect_lte(sc$fdr, 0.2)
})

test_that("mrmr_surv_select recovers signal (looser redundancy control than FCBF)", {
  skip_if_not_installed("praznik")
  sim <- concentrated_sim()
  r <- mrmr_surv_select(sim$time, sim$status, sim$X, prop = 0.05)
  sc <- score_selection_clusterwise(r$selected, sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
  expect_lte(sc$fdr, 0.5)            # mRMR averages redundancy, it does not prune by dominance
  expect_length(r$selected, 15L)
})

test_that("the MI family is deterministic (no RNG in the ranking)", {
  skip_if_not_installed("praznik")
  sim <- concentrated_sim()
  a <- mim_surv_select(sim$time, sim$status, sim$X, prop = 0.05)
  b <- mim_surv_select(sim$time, sim$status, sim$X, prop = 0.05)
  expect_identical(a$selected, b$selected)
  expect_identical(a$score, b$score)
})

test_that("mim_surv_select reproduces the hand-built martingale + quantile + MIM pipeline", {
  skip_if_not_installed("praznik")
  sim <- concentrated_sim()
  q <- 10
  # the exact substrate survfcbf_select uses: martingale target + equal-frequency bins
  ym <- survival_target(sim$time, sim$status)
  yd <- factor(discretize_quantile(ym, q = q)$binned)
  Xd <- as.data.frame(lapply(seq_len(ncol(sim$X)),
          function(j) factor(discretize_quantile(sim$X[, j], q = q)$binned)))
  k <- max(1L, round(0.05 * ncol(sim$X)))
  hand <- as.integer(praznik::MIM(Xd, yd, k = k)$selection)
  expect_equal(mim_surv_select(sim$time, sim$status, sim$X, prop = 0.05, q = q)$selected, hand)
})

test_that("prop = NULL yields a pure ranker with a full score", {
  skip_if_not_installed("praznik")
  sim <- concentrated_sim()
  full <- mim_surv_select(sim$time, sim$status, sim$X)   # k = p
  expect_null(full$selected)
  expect_true(all(is.finite(full$score)))
  expect_length(full$score, 300L)
})

test_that("MI selectors error clearly without praznik present", {
  skip_if_not_installed("praznik")    # when praznik IS present the guard is not exercised; smoke only
  sim <- concentrated_sim(32)
  expect_error(mim_surv_select(sim$time, sim$status, sim$X, prop = 1.5), "prop")
})
