# Acceptance tests for make_censoring() / calibrate_censoring().

sim_censor <- function(n = 10000, seed = 1) {
  bl <- make_sigma_blocks(p = 60, k_drivers = 4, cluster_size = 3, rho = 0.8,
                          n_distractor_subblocks = 3, distractor_size = 4,
                          distractor_rho = 0.7, interaction = TRUE)
  X <- sample_features(n, bl$Sigma, marginals = "gaussian", seed = seed, return_cor = FALSE)$X
  eta <- build_linear_predictor(X, bl$manifest, beta_scale = 1.2)
  Tt <- bender_weibull_time(eta, seed = seed + 100L)
  list(bl = bl, X = X, eta = eta, Tt = Tt)
}

test_that("realised censoring fraction matches the target across mechanisms and levels", {
  s <- sim_censor()
  for (mech in c("independent", "informative_eta", "informative_disjoint")) {
    for (f in c(0.2, 0.5, 0.7)) {
      param <- calibrate_censoring(s$eta, s$X, s$Tt, target_fraction = f, mechanism = mech,
                                   manifest = s$bl$manifest, seed = 7)
      C <- make_censoring(s$eta, s$X, param, mechanism = mech, manifest = s$bl$manifest, seed = 8)
      expect_true(abs(mean(C < s$Tt) - f) < 0.02, info = paste(mech, f))
    }
  }
})

test_that("administrative censoring hits the target fraction via the time quantile", {
  s <- sim_censor()
  for (f in c(0.2, 0.5, 0.7)) {
    tau <- calibrate_censoring(s$eta, s$X, s$Tt, target_fraction = f,
                               mechanism = "independent", fraction_model = "administrative")
    C <- make_censoring(s$eta, s$X, tau, mechanism = "independent",
                        fraction_model = "administrative")
    expect_true(abs(mean(C < s$Tt) - f) < 0.01, info = paste("administrative", f))
  }
})

test_that("independent censoring is unassociated with eta and covariates; informative arms associate", {
  s <- sim_censor()
  dcols <- default_disjoint_cols(s$bl$manifest)
  z <- rowMeans(s$X[, dcols])

  p_ind <- calibrate_censoring(s$eta, s$X, s$Tt, 0.5, "independent", manifest = s$bl$manifest, seed = 9)
  C_ind <- make_censoring(s$eta, s$X, p_ind, "independent", manifest = s$bl$manifest, seed = 10)
  expect_lt(abs(stats::cor(C_ind, s$eta)), 0.05)
  expect_lt(abs(stats::cor(C_ind, z)), 0.05)

  p_eta <- calibrate_censoring(s$eta, s$X, s$Tt, 0.5, "informative_eta", manifest = s$bl$manifest, seed = 9)
  C_eta <- make_censoring(s$eta, s$X, p_eta, "informative_eta", manifest = s$bl$manifest, seed = 10)
  expect_lt(stats::cor(C_eta, s$eta), -0.1)                 # higher eta -> shorter censoring time

  p_dis <- calibrate_censoring(s$eta, s$X, s$Tt, 0.5, "informative_disjoint", manifest = s$bl$manifest, seed = 9)
  C_dis <- make_censoring(s$eta, s$X, p_dis, "informative_disjoint", manifest = s$bl$manifest, seed = 10)
  expect_lt(stats::cor(C_dis, z), -0.1)                     # associates with the disjoint covariate
  expect_lt(abs(stats::cor(C_dis, s$eta)), 0.05)            # but not with eta
})

test_that("censoring draws are reproducible under a fixed seed", {
  s <- sim_censor(n = 1000, seed = 2)
  a <- make_censoring(s$eta, s$X, 0.5, "informative_eta", manifest = s$bl$manifest, seed = 11)
  b <- make_censoring(s$eta, s$X, 0.5, "informative_eta", manifest = s$bl$manifest, seed = 11)
  expect_identical(a, b)
})

test_that("invalid arguments are rejected", {
  s <- sim_censor(n = 500, seed = 3)
  expect_error(calibrate_censoring(s$eta, s$X, s$Tt, target_fraction = 0), "in \\(0, 1\\)")
  expect_error(calibrate_censoring(s$eta, s$X, s$Tt, 0.5, "informative_eta",
                                   fraction_model = "administrative"), "only mechanism")
  X0 <- s$X
  # informative_disjoint without manifest or disjoint_cols cannot find the subblock
  expect_error(make_censoring(s$eta, X0, 0.5, "informative_disjoint"), "manifest")
})
