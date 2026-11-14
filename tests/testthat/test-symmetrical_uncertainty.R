# Acceptance tests for the entropy / SU primitives.

test_that("entropy matches known values", {
  expect_equal(entropy(c(0, 0, 1, 1)), 1)        # uniform binary = 1 bit
  expect_equal(entropy(rep(1, 10)), 0)           # constant = 0
  expect_equal(entropy(c(1, 2, 3, 4)), 2)        # uniform over 4 categories = 2 bits
})

test_that("conditional entropy is zero for self and bounded by the marginal", {
  x <- c(1, 1, 2, 2, 3, 3)
  y <- c(1, 2, 1, 2, 1, 2)
  expect_equal(conditional_entropy(x, x), 0)                       # H(X | X) = 0
  expect_lte(conditional_entropy(x, y), entropy(x) + 1e-12)        # H(X | Y) <= H(X)
})

test_that("symmetrical uncertainty is in [0, 1], symmetric, 1 for self, 0 for constant", {
  x <- c(1, 1, 2, 2, 1, 2, 1, 2)
  y <- c(1, 2, 1, 2, 2, 1, 2, 1)
  expect_equal(symmetrical_uncertainty(x, x), 1)
  expect_equal(symmetrical_uncertainty(x, rep(0, 8)), 0)
  expect_equal(symmetrical_uncertainty(x, y), symmetrical_uncertainty(y, x))
  su <- symmetrical_uncertainty(x, y)
  expect_gte(su, 0)
  expect_lte(su, 1)
})

test_that("symmetrical uncertainty matches hand-computed values", {
  # independent: x = (0,0,1,1), y = (0,1,0,1) -> IG = 0 -> SU = 0
  expect_equal(symmetrical_uncertainty(c(0, 0, 1, 1), c(0, 1, 0, 1)), 0)

  # x = (0,0,0,1), w = (0,0,1,1):
  #   H(x) = 0.8113, H(w) = 1, H(x|w) = 0.5, IG = 0.3113
  #   SU = 2 * 0.3113 / (0.8113 + 1) = 0.3437
  expect_equal(symmetrical_uncertainty(c(0, 0, 0, 1), c(0, 0, 1, 1)), 0.3437, tolerance = 1e-3)
})

test_that("primitives reject NA and mismatched lengths", {
  expect_error(entropy(c(1, NA, 2)), "without NA")
  expect_error(conditional_entropy(c(1, 2, 3), c(1, 2)), "equal-length")
})

test_that("uniform weights reproduce the unweighted primitives", {
  set.seed(1)
  x <- sample(1:3, 100, replace = TRUE)
  y <- sample(1:2, 100, replace = TRUE)
  w <- rep(2, 100)
  expect_equal(entropy(x, w), entropy(x))
  expect_equal(conditional_entropy(x, y, w), conditional_entropy(x, y))
  expect_equal(symmetrical_uncertainty(x, y, w), symmetrical_uncertainty(x, y))
})

test_that("integer weights equal observation replication", {
  x <- c(1, 1, 2, 2, 2)
  y <- c(1, 2, 1, 2, 2)
  w <- c(2, 2, 1, 1, 1)
  idx <- rep(seq_along(x), times = w)            # replicate observations by their weight
  expect_equal(entropy(x, w), entropy(x[idx]))
  expect_equal(symmetrical_uncertainty(x, y, w), symmetrical_uncertainty(x[idx], y[idx]))
})

test_that("weighted primitives reject malformed weights", {
  expect_error(entropy(c(1, 2, 3), w = c(1, 2)), "length")
  expect_error(entropy(c(1, 2, 3), w = c(-1, 1, 1)), "non-negative")
})
