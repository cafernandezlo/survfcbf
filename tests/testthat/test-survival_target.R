# Acceptance tests for survival_target().

test_that("martingale residuals sum to ~0 and have the right length", {
  set.seed(1)
  n <- 200
  time <- stats::rexp(n, rate = 1)
  status <- rbinom(n, 1, 0.6)
  ym <- survival_target(time, status)
  expect_length(ym, n)
  expect_equal(sum(ym), 0, tolerance = 1e-6)        # martingale residuals sum to zero
})

test_that("an event gives a higher residual than a censored subject at the same time", {
  # two subjects share each time; the event one should sit exactly 1 above the censored one
  time <- c(2, 2, 5, 5, 9, 9)
  status <- c(1, 0, 1, 0, 1, 0)
  ym <- survival_target(time, status)
  event <- ym[status == 1]
  cens <- ym[status == 0]
  expect_true(all(event > cens))                    # status raises the residual
  expect_equal(event - cens, rep(1, 3), tolerance = 1e-8)
})

test_that("early events get more positive residuals than long survivors", {
  set.seed(2)
  n <- 300
  time <- stats::rexp(n, 1)
  status <- rep(1L, n)                               # all events
  ym <- survival_target(time, status)
  # residual decreases with time (early event -> high, late event -> low/negative)
  expect_lt(stats::cor(time, ym), -0.8)
})

test_that("survival_target is deterministic and validates its inputs", {
  time <- c(1, 2, 3, 4); status <- c(1, 0, 1, 1)
  expect_identical(survival_target(time, status), survival_target(time, status))
  expect_identical(survival_target(time, as.logical(status)), survival_target(time, status))

  expect_error(survival_target(c(1, -1, 2), c(1, 0, 1)), "positive")
  expect_error(survival_target(c(1, 2, 3), c(1, 2, 0)), "0/1")
  expect_error(survival_target(c(1, 2, 3), c(1, 0)), "same length")
})
