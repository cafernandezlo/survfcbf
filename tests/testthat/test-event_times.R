# Acceptance tests for the event-time arms.

test_that("Weibull PH: baseline median is 1 at eta = 0", {
  Tt <- bender_weibull_time(rep(0, 50000), shape_a = 1.5, scale_lambda = log(2), seed = 1)
  expect_equal(stats::median(Tt), 1, tolerance = 0.02)
})

test_that("Weibull PH: higher eta gives stochastically shorter times", {
  T0 <- bender_weibull_time(rep(0, 20000), seed = 2)
  T1 <- bender_weibull_time(rep(1, 20000), seed = 3)
  expect_lt(stats::median(T1), stats::median(T0))
  expect_lt(mean(T1 > 1), mean(T0 > 1))            # less mass above the baseline median
})

test_that("Weibull PH: a = 1 reduces to the exponential", {
  Tt <- bender_weibull_time(rep(0, 20000), shape_a = 1, scale_lambda = log(2), seed = 4)
  ks <- suppressWarnings(stats::ks.test(Tt, "pexp", rate = log(2))$statistic)
  expect_lt(as.numeric(ks), 0.03)
})

test_that("Weibull PH: draws are reproducible under a fixed seed", {
  a <- bender_weibull_time(rep(0.3, 1000), seed = 5)
  b <- bender_weibull_time(rep(0.3, 1000), seed = 5)
  expect_identical(a, b)
})

test_that("log-logistic AFT: baseline is log-logistic with median 1 at eta = 0", {
  Tt <- aft_loglogistic_time(rep(0, 60000), shape = 3, seed = 6)
  expect_equal(stats::median(Tt), 1, tolerance = 0.02)
  # log-logistic quantiles: q(p) = (p/(1-p))^(1/shape)
  expect_equal(as.numeric(stats::quantile(Tt, 0.75)), (0.75 / 0.25)^(1 / 3), tolerance = 0.03)
  expect_equal(as.numeric(stats::quantile(Tt, 0.25)), (0.25 / 0.75)^(1 / 3), tolerance = 0.03)
})

test_that("log-logistic AFT: higher eta accelerates failure (AFT time-scaling)", {
  T0 <- aft_loglogistic_time(rep(0, 40000), shape = 3, seed = 7)
  T1 <- aft_loglogistic_time(rep(1, 40000), shape = 3, seed = 8)
  expect_lt(stats::median(T1), stats::median(T0))
  # median ratio reflects exp(-delta eta) time scaling
  expect_equal(stats::median(T1) / stats::median(T0), exp(-1), tolerance = 0.05)
})

test_that("the AFT arm violates proportional hazards while the Weibull arm respects it", {
  skip_if_not_installed("survival")
  n <- 4000
  g <- rep(c(0, 1), each = n)
  eta <- 1.2 * g
  status <- rep(1L, length(g))

  Tw <- bender_weibull_time(eta, shape_a = 1.5, scale_lambda = log(2), seed = 20)
  pw <- survival::cox.zph(survival::coxph(survival::Surv(Tw, status) ~ g))$table["g", "p"]

  Ta <- aft_loglogistic_time(eta, shape = 3, seed = 21)
  pa <- survival::cox.zph(survival::coxph(survival::Surv(Ta, status) ~ g))$table["g", "p"]

  expect_lt(pa, 0.01)        # log-logistic AFT strongly violates PH
  expect_gt(pw, 0.05)        # Weibull PH does not
})

test_that("event-time arms reject malformed arguments", {
  expect_error(bender_weibull_time(c(0, 1), shape_a = -1), "positive")
  expect_error(bender_weibull_time(c(0, NA)), "finite numeric")
  expect_error(aft_loglogistic_time(c(0, 1), shape = 0), "positive")
})
