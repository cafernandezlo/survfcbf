# Acceptance tests for calibrate_beta().

blocks24 <- function() {
  make_sigma_blocks(p = 24, k_drivers = 3, cluster_size = 3, rho = 0.8,
                    n_distractor_subblocks = 2, distractor_size = 4,
                    distractor_rho = 0.7, interaction = TRUE)
}

# analytic pairwise concordance of a linear predictor (what calibrate_beta targets)
pair_C <- function(eta, n_pairs = 3e5, seed = 99) {
  set.seed(seed)
  n <- length(eta)
  i <- sample.int(n, n_pairs, replace = TRUE)
  j <- sample.int(n, n_pairs, replace = TRUE)
  same <- i == j
  j[same] <- (i[same] %% n) + 1L
  mean(stats::plogis(abs(eta[i] - eta[j])))
}

test_that("calibrated beta_scale hits the target analytic concordance within tol", {
  bl <- blocks24()
  target <- 0.70
  bs <- calibrate_beta(bl$Sigma, bl$manifest, target_C = target,
                       marginals = "gaussian", tol = 0.003, seed = 1)
  X <- sample_features(20000, bl$Sigma, marginals = "gaussian", seed = 2, return_cor = FALSE)$X
  eta <- build_linear_predictor(X, bl$manifest, beta_scale = bs)
  expect_equal(pair_C(eta), target, tolerance = 0.01)
})

test_that("calibrated scale matches empirical Harrell's C on generated PH times", {
  skip_if_not_installed("survival")
  bl <- blocks24()
  target <- 0.70
  bs <- calibrate_beta(bl$Sigma, bl$manifest, target_C = target,
                       marginals = "gaussian", tol = 0.003, seed = 3)
  X <- sample_features(20000, bl$Sigma, marginals = "gaussian", seed = 4, return_cor = FALSE)$X
  eta <- build_linear_predictor(X, bl$manifest, beta_scale = bs)
  set.seed(5)
  Tt <- stats::rexp(length(eta), rate = exp(eta))      # PH exponential times, uncensored
  # survival::concordance treats a larger RHS as predicting longer survival; eta is
  # a risk score (larger -> earlier event), so Harrell's risk C-index is 1 - cidx.
  cidx <- survival::concordance(survival::Surv(Tt, rep(1L, length(eta))) ~ eta)$concordance
  expect_equal(1 - cidx, target, tolerance = 0.02)
})

test_that("concordance is monotone in beta_scale (higher target -> larger scale)", {
  bl <- blocks24()
  bs_lo <- calibrate_beta(bl$Sigma, bl$manifest, target_C = 0.60, marginals = "gaussian", seed = 6)
  bs_hi <- calibrate_beta(bl$Sigma, bl$manifest, target_C = 0.80, marginals = "gaussian", seed = 6)
  expect_gt(bs_hi, bs_lo)
})

test_that("calibration is reproducible under a fixed seed", {
  bl <- blocks24()
  a <- calibrate_beta(bl$Sigma, bl$manifest, target_C = 0.7, marginals = "gaussian", seed = 8)
  b <- calibrate_beta(bl$Sigma, bl$manifest, target_C = 0.7, marginals = "gaussian", seed = 8)
  expect_equal(a, b)
})

test_that("invalid target_C is rejected", {
  bl <- blocks24()
  expect_error(calibrate_beta(bl$Sigma, bl$manifest, target_C = 0.4), "0.5, 1")
  expect_error(calibrate_beta(bl$Sigma, bl$manifest, target_C = 1.2), "0.5, 1")
})
