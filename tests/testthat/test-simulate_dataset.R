# Acceptance tests for simulate_dataset().
# Uses a reduced p for speed; the event count and diagnostics are p-independent.

test_that("simulate_dataset composes the pipeline and reports coherent diagnostics", {
  pr <- default_params(p = 300, n = 1000, censor_fraction = 0.5, target_C = 0.70,
                       rho = 0.8, marginals = "gaussian", seed = 1, calib_n_presample = 2000)
  sim <- simulate_dataset(pr)

  expect_equal(dim(sim$X), c(1000L, 300L))
  expect_length(sim$time, 1000L)
  expect_true(all(sim$status %in% c(0L, 1L)))
  expect_true(all(sim$time > 0))
  expect_equal(nrow(sim$manifest), 300L)
  expect_setequal(names(sim$realized),
                  c("beta_scale", "censoring_fraction", "n_events", "concordance", "intracluster_cor"))

  r <- sim$realized
  expect_true(abs(r$censoring_fraction - 0.5) < 0.05)        # calibrated censoring fraction
  expect_equal(r$concordance, 0.70, tolerance = 0.03)        # PH arm matches the calibrated target
  expect_equal(r$intracluster_cor, 0.8, tolerance = 0.05)    # gaussian -> latent rho
  expect_equal(r$n_events, sum(sim$status))
})

test_that("the default n / censoring yields >= 50 effective events", {
  sim <- simulate_dataset(default_params(p = 300, n = 200, censor_fraction = 0.5, seed = 2,
                                         calib_n_presample = 1500))
  expect_gt(sim$realized$n_events, 50)
})

test_that("simulation is reproducible under a fixed seed and varies otherwise", {
  pr <- default_params(p = 300, n = 150, seed = 42, calib_n_presample = 1200)
  a <- simulate_dataset(pr)
  b <- simulate_dataset(pr)
  expect_identical(a$X, b$X)
  expect_identical(a$time, b$time)
  expect_identical(a$status, b$status)

  c2 <- simulate_dataset(default_params(p = 300, n = 150, seed = 99, calib_n_presample = 1200))
  expect_false(identical(a$time, c2$time))
})

test_that("the non-PH / informative arm runs end to end with sane output", {
  pr <- default_params(p = 300, n = 800, time_model = "loglogistic_aft",
                       censor_mechanism = "informative_eta", censor_fraction = 0.6,
                       seed = 3, calib_n_presample = 1500)
  sim <- simulate_dataset(pr)

  expect_equal(dim(sim$X), c(800L, 300L))
  expect_true(all(sim$status %in% c(0L, 1L)))
  expect_true(abs(sim$realized$censoring_fraction - 0.6) < 0.05)
  expect_equal(sim$realized$concordance, 0.70, tolerance = 0.04)   # AFT arm calibrated to target_C (empirical bisection)
})

test_that("the AFT arm tracks target_C rather than the PH-formula overshoot", {
  # discriminating regression: the analytic plogis calibration overshot here
  # (target 0.65 realised ~0.80); the empirical AFT calibration tracks the target.
  cs <- vapply(c(11, 12, 13), function(s) {
    pr <- default_params(p = 300, n = 700, time_model = "loglogistic_aft",
                         target_C = 0.65, censor_fraction = 0.3, seed = s,
                         calib_n_presample = 1500)
    simulate_dataset(pr)$realized$concordance
  }, numeric(1))
  expect_equal(mean(cs), 0.65, tolerance = 0.03)
  expect_true(all(cs < 0.72))                                # well below the pre-fix ~0.80 overshoot
})

test_that("default_params applies overrides and leaves the rest at ficha defaults", {
  pr <- default_params(n = 400L, censor_fraction = 0.7)
  expect_equal(pr$n, 400L)
  expect_equal(pr$censor_fraction, 0.7)
  expect_equal(pr$p, 5000L)             # untouched default
  expect_equal(pr$target_C, 0.70)
})

test_that("low-variance-signal scenario shrinks driver-cluster variance but preserves structure", {
  base <- list(p = 300, n = 1000, seed = 7, calib_n_presample = 1500)
  hi <- simulate_dataset(do.call(default_params, c(base, list(driver_var_scale = 1))))
  lo <- simulate_dataset(do.call(default_params, c(base, list(driver_var_scale = 0.2))))

  m <- hi$manifest
  sig <- which(m$role %in% c("driver", "redundant"))
  noise <- which(m$role == "noise")
  v_lo <- mean(apply(lo$X[, sig], 2, stats::var))
  v_hi <- mean(apply(hi$X[, sig], 2, stats::var))

  expect_lt(v_lo, 0.1 * v_hi)                                   # driver-cluster variance collapses
  expect_lt(v_lo, mean(apply(lo$X[, noise], 2, stats::var)))    # below noise: a variance filter is misled

  # rank structure and outcome are untouched (scaling is post-eta and scale-invariant for cor)
  expect_equal(lo$realized$intracluster_cor, hi$realized$intracluster_cor, tolerance = 1e-8)
  expect_identical(lo$status, hi$status)
})
