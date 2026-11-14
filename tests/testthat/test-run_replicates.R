# Acceptance tests for run_replicates().

test_that("run_replicates summarises diagnostics with Monte Carlo standard errors", {
  pr <- default_params(p = 300, n = 200, censor_fraction = 0.5, target_C = 0.70,
                       calib_n_presample = 1500)
  out <- run_replicates(pr, R = 12, master_seed = 1)

  expect_equal(out$R, 12L)
  expect_setequal(out$summary$diagnostic,
                  c("beta_scale", "censoring_fraction", "n_events", "concordance", "intracluster_cor"))
  expect_true(all(out$summary$mc_se >= 0))
  expect_equal(dim(out$replicates), c(5L, 12L))

  pick <- function(k) out$summary$mean[out$summary$diagnostic == k]
  expect_true(abs(pick("censoring_fraction") - 0.5) < 0.05)   # averaged near the target fraction
  expect_true(abs(pick("concordance") - 0.70) < 0.03)         # averaged near target_C
  expect_true(abs(pick("intracluster_cor") - 0.8) < 0.04)     # gaussian -> latent rho

  # mc_se equals sd / sqrt(R)
  se_manual <- stats::sd(out$replicates["concordance", ]) / sqrt(12)
  expect_equal(out$summary$mc_se[out$summary$diagnostic == "concordance"], se_manual, tolerance = 1e-10)
})

test_that("run_replicates is reproducible from the master seed and varies with it", {
  pr <- default_params(p = 300, n = 150, calib_n_presample = 1000)
  a <- run_replicates(pr, R = 6, master_seed = 42)
  b <- run_replicates(pr, R = 6, master_seed = 42)
  expect_identical(a$replicates, b$replicates)

  c2 <- run_replicates(pr, R = 6, master_seed = 99)
  expect_false(identical(a$replicates, c2$replicates))
})
