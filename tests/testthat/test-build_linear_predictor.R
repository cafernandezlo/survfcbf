# Acceptance tests for build_linear_predictor().

sim_eta <- function(beta_scale = 1, beta_int = 1, n = 4000, seed = 1) {
  bl <- make_sigma_blocks(p = 24, k_drivers = 3, cluster_size = 3, rho = 0.8,
                          n_distractor_subblocks = 2, distractor_size = 4,
                          distractor_rho = 0.7, interaction = TRUE)
  sf <- sample_features(n = n, Sigma = bl$Sigma, marginals = "gaussian",
                        seed = seed, return_cor = FALSE)
  list(X = sf$X, manifest = bl$manifest,
       eta = build_linear_predictor(sf$X, bl$manifest, beta_scale, beta_int))
}

test_that("eta is constant when beta_scale = 0", {
  s <- sim_eta(beta_scale = 0)
  expect_equal(stats::var(s$eta), 0)
})

test_that("eta has one entry per observation", {
  s <- sim_eta(n = 500)
  expect_length(s$eta, 500L)
})

test_that("interaction features have ~0 marginal but nonzero joint association with eta", {
  s <- sim_eta(beta_scale = 1, beta_int = 1, n = 6000, seed = 2)
  inter <- which(s$manifest$role == "interaction")
  Xa <- s$X[, inter[1]]; Xb <- s$X[, inter[2]]

  expect_lt(abs(stats::cor(Xa, s$eta)), 0.08)          # marginal ~ 0
  expect_lt(abs(stats::cor(Xb, s$eta)), 0.08)

  prod_term <- as.numeric(scale(Xa)) * as.numeric(scale(Xb))
  expect_gt(abs(stats::cor(prod_term, s$eta)), 0.2)    # joint nonzero
})

test_that("with beta_int = 0 eta is exactly the additive driver effect", {
  s <- sim_eta(beta_scale = 1, beta_int = 0, n = 3000, seed = 3)
  signal <- which(s$manifest$beta != 0)
  recon <- as.vector(s$X[, signal] %*% s$manifest$beta[signal])
  expect_equal(s$eta, recon)
})

test_that("beta_scale linearly scales the whole linear predictor", {
  s1 <- sim_eta(beta_scale = 1, beta_int = 1, n = 2000, seed = 4)
  s2 <- sim_eta(beta_scale = 2, beta_int = 1, n = 2000, seed = 4)
  expect_equal(s2$eta, 2 * s1$eta)
})

test_that("input guards reject malformed arguments", {
  bl <- make_sigma_blocks(p = 24, k_drivers = 3, cluster_size = 3, rho = 0.8,
                          n_distractor_subblocks = 2, distractor_size = 4,
                          distractor_rho = 0.7)
  X <- matrix(0, nrow = 5, ncol = 24)
  expect_error(build_linear_predictor(X, bl$manifest[1:10, ], beta_scale = 1),
               "one row per column")
  expect_error(build_linear_predictor(X, bl$manifest, beta_scale = c(1, 2)),
               "single finite number")
})
