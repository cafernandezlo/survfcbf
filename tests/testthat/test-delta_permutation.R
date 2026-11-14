# Acceptance tests for delta_permutation().

su_relevance <- function(X, y) {
  vapply(seq_len(ncol(X)), function(j) symmetrical_uncertainty(X[, j], y), numeric(1))
}

test_that("under a null target few/no features exceed the permutation delta", {
  set.seed(5)
  n <- 150; p <- 40
  y <- sample(1:3, n, replace = TRUE)
  X <- matrix(sample(1:3, n * p, replace = TRUE), n, p)   # all independent of y

  delta <- delta_permutation(X, y, B = 100, prob = 0.95, seed = 1)
  expect_true(delta > 0 && delta < 1)
  n_exceed <- sum(su_relevance(X, y) >= delta)
  expect_lte(n_exceed, 3L)                                  # ~ family-wise control near 5%
})

test_that("a feature strongly associated with the target exceeds delta", {
  set.seed(6)
  n <- 150; p <- 20
  y <- sample(1:3, n, replace = TRUE)
  X <- matrix(sample(1:3, n * p, replace = TRUE), n, p)
  X[, 1] <- y                                               # planted perfect signal

  delta <- delta_permutation(X, y, B = 100, prob = 0.95, seed = 2)
  su <- su_relevance(X, y)
  expect_equal(su[1], 1)                                    # SU(y, y) = 1
  expect_gt(su[1], delta)                                   # the signal is selectable
})

test_that("delta is reproducible under a seed and increases with prob", {
  set.seed(7)
  n <- 120; p <- 25
  y <- sample(1:2, n, replace = TRUE)
  X <- matrix(sample(1:2, n * p, replace = TRUE), n, p)

  a <- delta_permutation(X, y, B = 80, prob = 0.95, seed = 3)
  b <- delta_permutation(X, y, B = 80, prob = 0.95, seed = 3)
  expect_equal(a, b)

  lo <- delta_permutation(X, y, B = 80, prob = 0.90, seed = 3)
  hi <- delta_permutation(X, y, B = 80, prob = 0.99, seed = 3)
  expect_lte(lo, hi)
})

test_that("delta_permutation validates its inputs", {
  X <- matrix(sample(1:2, 60, replace = TRUE), 20, 3)
  expect_error(delta_permutation(X, rep(1, 19)), "length nrow")
  expect_error(delta_permutation(X, rep(1, 20), prob = 1), "in \\(0, 1\\)")
})
