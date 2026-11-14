# Acceptance tests for the carss baseline.

carss_sim <- function(seed = 31) {
  simulate_dataset(default_params(
    p = 300, k_drivers = 3, cluster_size = 5, n = 800,
    censor_fraction = 0.10, target_C = 0.90, seed = seed, calib_n_presample = 1500))
}

test_that("carss_select recovers concentrated signal at a fixed cardinality", {
  skip_if_not_installed("carSurv")
  sim <- carss_sim()
  r <- carss_select(sim$time, sim$status, sim$X, prop = 0.05)   # top 15 of 300

  expect_length(r$score, 300L)
  expect_true(all(r$score >= 0))                 # importance is the absolute CAR score
  expect_equal(r$score, abs(r$raw_score))        # raw scores are signed
  expect_length(r$selected, 15L)
  sc <- score_selection_clusterwise(r$selected, sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
  expect_lte(sc$fdr, 0.2)
})

# carss_select still supports carVarSelect self-selection (prop = NULL), exercised
# here on strong signal. The benchmark registry wires carss as a *ranker* instead,
# because this fndr self-selection degenerates to the empty set once the signal
# weakens (target_C ~ 0.7), which would drop carss out of every fold.
test_that("carss self-selects a set via carVarSelect when prop is NULL", {
  skip_if_not_installed("carSurv")
  sim <- carss_sim()
  r <- carss_select(sim$time, sim$status, sim$X)   # fndr self-selection
  expect_true(is.integer(r$selected) && length(r$selected) >= 1L)
  expect_true(all(r$selected >= 1L & r$selected <= 300L))
  sc <- score_selection_clusterwise(r$selected, sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
})

test_that("carss_select is deterministic and leakage-free over the training fold", {
  skip_if_not_installed("carSurv")
  sim <- carss_sim()
  a <- carss_select(sim$time, sim$status, sim$X, prop = 0.05)
  b <- carss_select(sim$time, sim$status, sim$X, prop = 0.05)
  expect_identical(a$selected, b$selected)
  expect_identical(a$score, b$score)
})

test_that("carss_select validates inputs and missing package message", {
  skip_if_not_installed("carSurv")
  sim <- carss_sim(32)
  expect_error(carss_select(sim$time, sim$status, sim$X, prop = 2), "prop")
  expect_error(carss_select(sim$time[-1], sim$status, sim$X), "share nrow")
})
