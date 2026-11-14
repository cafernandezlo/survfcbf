test_that("bayes_signed_rank detects a clear win, a tie and the ROPE", {
  set.seed(1)
  a <- c(0.62, 0.58, 0.66, 0.61, 0.59, 0.64, 0.60)
  # x clearly better: every pair differs by well over the ROPE
  b <- a - 0.10
  r <- bayes_signed_rank(a, b, rope = 0.02, nsim = 2000)
  expect_equal(unname(r["n"]), 7)
  expect_gt(r[["p_x"]], 0.9)
  expect_lt(r[["p_y"]], 0.05)

  # identical scores: all the mass is practical equivalence
  r0 <- bayes_signed_rank(a, a, rope = 0.02, nsim = 2000)
  expect_gt(r0[["p_rope"]], 0.99)

  # probabilities form a distribution
  expect_equal(sum(r[c("p_x", "p_rope", "p_y")]), 1, tolerance = 1e-8)
})

test_that("argmax is the default summary", {
  set.seed(3)
  a <- c(0.62, 0.58, 0.66, 0.61, 0.59, 0.64, 0.60)
  b <- a - c(0.03, 0.02, 0.04, 0.03, 0.01, 0.03, 0.02)
  set.seed(4); d <- bayes_signed_rank(a, b, rope = 0.01, nsim = 1000)
  set.seed(4); g <- bayes_signed_rank(a, b, rope = 0.01, nsim = 1000, summary = "argmax")
  expect_identical(d, g)
})

test_that("the two summaries of the posterior differ, and only 'mean' has a ceiling", {
  # The "mean" summary is bounded by P((Z1+Z2)/2 > rope); Benavoli's "argmax"
  # convention (the default) is not, and converges to 1. A small real effect separates
  # them sharply.
  set.seed(11)
  z <- rnorm(120, mean = 0.03, sd = 0.055)
  zero <- rep(0, length(z))
  ceiling_mean <- 1 - pnorm(0.02, 0.03, 0.055 / sqrt(2))          # ~0.60

  m <- bayes_signed_rank(z, zero, rope = 0.02, nsim = 3000, summary = "mean")
  a <- bayes_signed_rank(z, zero, rope = 0.02, nsim = 3000, summary = "argmax")

  expect_lt(m[["p_x"]], ceiling_mean + 0.1)                      # stuck near the ceiling
  expect_gt(a[["p_x"]], 0.95)                                    # decisive under argmax
  expect_equal(sum(a[c("p_x", "p_rope", "p_y")]), 1, tolerance = 1e-8)
  expect_equal(unname(m["n"]), unname(a["n"]))                    # same data, same pairs
})

test_that("bayes_signed_rank honours higher_better and drops incomplete pairs", {
  set.seed(1)
  a <- c(0.20, 0.18, 0.25, 0.22, 0.19)
  b <- a + 0.10                                   # b has the larger (worse, for a loss) values
  hi <- bayes_signed_rank(a, b, rope = 0.02, higher_better = TRUE,  nsim = 2000)
  lo <- bayes_signed_rank(a, b, rope = 0.02, higher_better = FALSE, nsim = 2000)
  expect_gt(hi[["p_y"]], 0.9)                     # larger is better -> b wins
  expect_gt(lo[["p_x"]], 0.9)                     # smaller is better -> a wins

  # NA pairs are dropped, not treated as ties
  aNA <- c(a, NA_real_); bNA <- c(b, 0.5)
  expect_equal(unname(bayes_signed_rank(aNA, bNA, rope = 0.02, nsim = 500)["n"]), 5)
})

test_that("bayes_signed_rank validates its inputs and needs two pairs", {
  expect_error(bayes_signed_rank(1:5, 1:4, rope = 0.02), "same length")
  expect_error(bayes_signed_rank(c(1, 2), c(1, 2), rope = 0), "positive number")
  r <- bayes_signed_rank(c(1, NA), c(2, 3), rope = 0.02, nsim = 100)
  expect_true(is.na(r[["p_x"]]))
})
