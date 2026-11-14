# Acceptance tests for discretize_quantile().

test_that("learned bins are approximately equal-frequency", {
  set.seed(1)
  x <- stats::rnorm(1000)
  d <- discretize_quantile(x, q = 5)
  counts <- tabulate(d$binned, nbins = 5)
  expect_length(counts, 5L)
  expect_true(all(counts >= 180 & counts <= 220))   # ~200 each, +/- 10%
})

test_that("re-applying training breaks uses no test information", {
  set.seed(2)
  train <- stats::rnorm(500)
  test <- stats::rnorm(300)
  learned <- discretize_quantile(train, q = 4)
  reapplied <- discretize_quantile(test, breaks = learned$breaks)
  expect_identical(reapplied$breaks, learned$breaks)        # same cutpoints, learned on train
  # binning test with train breaks matches a manual cut with those breaks
  expect_equal(reapplied$binned,
               as.integer(cut(test, learned$breaks, labels = FALSE, include.lowest = TRUE)))
})

test_that("out-of-range test values fall in the extreme bins", {
  train <- seq(0, 1, length.out = 200)
  learned <- discretize_quantile(train, q = 5)
  edge <- discretize_quantile(c(-100, 100), breaks = learned$breaks)$binned
  expect_equal(edge[1], 1L)                  # below the training range -> first bin
  expect_equal(edge[2], 5L)                  # above the training range -> last bin
})

test_that("the default q follows Bommert's rule", {
  set.seed(3)
  x <- stats::rnorm(200)                      # floor(200/3)=66 -> min(66,10)=10 -> q=10
  d <- discretize_quantile(x)
  expect_equal(length(unique(d$binned)), 10L)
  expect_equal(max(d$binned), 10L)
})

test_that("few distinct values collapse to fewer bins without error", {
  x <- rep(c(1, 2, 3), times = c(40, 30, 30))
  d <- discretize_quantile(x, q = 10)
  expect_lte(length(unique(d$binned)), 3L)    # cannot exceed the number of distinct values
  expect_false(anyNA(d$binned))
})

test_that("discretization is deterministic and rejects NA", {
  x <- stats::qnorm(stats::ppoints(50))
  expect_identical(discretize_quantile(x, q = 4), discretize_quantile(x, q = 4))
  expect_error(discretize_quantile(c(1, NA, 2)), "without NA")
})
