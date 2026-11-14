# Acceptance tests for fcbf_select().

test_that("a tiebreak key orders features tied on relevance, and top-K stays an exact prefix", {
  set.seed(5)
  n <- 120
  y <- sample.int(2, n, TRUE)
  base <- ifelse(runif(n) < 0.3, sample.int(2, n, TRUE), y)
  X <- cbind(a = base, b = base, c = sample.int(2, n, TRUE), d = sample.int(2, n, TRUE))
  tb <- c(0.1, 0.9, 0.5, 0.2)                                 # a and b tie on relevance
  expect_identical(colnames(X)[fcbf_select(X, y)$selected][1], "a")                 # column order
  expect_identical(colnames(X)[fcbf_select(X, y, tiebreak = tb)$selected][1], "b")  # the key
  full <- fcbf_select(X, y, tiebreak = tb)$selected
  for (K in seq_along(full))
    expect_identical(fcbf_select(X, y, tiebreak = tb, top = K)$selected, utils::head(full, K))
  expect_error(fcbf_select(X, y, tiebreak = 1:3), "tiebreak")
  expect_error(fcbf_select(X, y, tiebreak = c(NA, 1, 1, 1)), "tiebreak")
})

test_that("fcbf_select keeps predominant features, prunes redundant ones, and traces decisions", {
  y <- c(0, 0, 0, 0, 1, 1, 1, 1)
  A <- c(0, 0, 1, 1, 1, 1, 1, 1)   # relevant
  B <- c(0, 0, 0, 0, 0, 0, 1, 1)   # relevant, not redundant with A
  C <- A                           # exact copy of A -> redundant, blanketed by A
  D <- c(0, 1, 0, 1, 0, 1, 0, 1)   # independent of y -> SU = 0
  X <- cbind(A = A, B = B, C = C, D = D)

  res <- fcbf_select(X, y, delta = 0.2)

  expect_equal(sort(res$selected), c(1L, 2L))                 # A and B are predominant
  expect_equal(res$trace$removed$feature, 3L)                 # C is pruned
  expect_equal(res$trace$removed$blanket, 1L)                 # blanketed by A
  expect_gte(res$trace$removed$su_blanket_feature, res$trace$removed$su_feature_target)
  expect_setequal(res$trace$representative_per_cluster[["1"]], c(1L, 3L))  # A represents {A, C}
  expect_setequal(res$trace$representative_per_cluster[["2"]], 2L)         # B represents {B}
  expect_equal(unname(res$trace$relevance_su[["D"]]), 0)      # D is irrelevant
})

test_that("fcbf_select uses precomputed relevance and makes y optional", {
  f1 <- rep(1:3, length.out = 120)
  X <- cbind(a = f1, b = f1, d = rep(1:2, length.out = 120))   # b is an exact copy of a
  res <- fcbf_select(X, delta = 0.05, relevance = c(0.9, 0.8, 0.2))   # no y supplied
  expect_equal(unname(res$trace$relevance_su), c(0.9, 0.8, 0.2))      # uses the precomputed relevance
  expect_true(1L %in% res$selected)
  expect_false(2L %in% res$selected)                                  # the redundant copy is pruned
  expect_true(2L %in% res$trace$removed$feature)
})

test_that("delta drops features below the relevance threshold", {
  y <- c(0, 0, 0, 0, 1, 1, 1, 1)
  A <- c(0, 0, 1, 1, 1, 1, 1, 1)   # SU(A, y) ~ 0.344
  D <- c(0, 1, 0, 1, 0, 1, 0, 1)   # SU = 0
  X <- cbind(A, D)

  expect_equal(fcbf_select(X, y, delta = 0.2)$selected, 1L)   # A passes, D filtered
  expect_length(fcbf_select(X, y, delta = 0.9)$selected, 0L)  # delta above SU(A, y): nothing selected
})

test_that("relevance ordering and self-redundancy: identical features collapse to one", {
  y <- c(0, 0, 1, 1, 0, 0, 1, 1)
  f <- c(0, 0, 1, 1, 1, 1, 1, 1)
  X <- cbind(f1 = f, f2 = f, f3 = f)        # three identical features
  res <- fcbf_select(X, y, delta = 0)

  expect_length(res$selected, 1L)                              # only one survives
  expect_equal(res$selected, 1L)                               # the first in rank order
  expect_setequal(res$trace$representative_per_cluster[["1"]], c(1L, 2L, 3L))
})

test_that("input guards", {
  X <- matrix(c(0, 1, 0, 1), nrow = 2)
  expect_error(fcbf_select(X, y = c(0, 1, 0), delta = 0), "length nrow")
  expect_error(fcbf_select(X, y = c(0, 1), delta = 1.5), "in \\[0, 1\\]")
  expect_error(fcbf_select(as.list(1:4), y = 1:2), "matrix or data.frame")
})

test_that("a top cap returns exactly the top-K predominant features of the full FCBF set", {
  # 4-class target with two orthogonal relevant contrasts (gA, gB), each with an
  # exact copy the Markov blanket must prune, plus a third correlated relevant one.
  y  <- rep(0:3, each = 50)
  gA <- as.integer(y >= 2)              # contrast 1
  gB <- as.integer(y %% 2 == 1)         # contrast 2, independent of gA -> neither blankets the other
  gC <- as.integer(y >= 1)              # contrast 3, relevant
  X  <- cbind(gA = gA, gB = gB, gA_copy = gA, gB_copy = gB, gC = gC)

  full <- fcbf_select(X, y, delta = 0)
  expect_gte(length(full$selected), 2L)                 # the orthogonal contrasts both survive
  expect_false(any(c(3L, 4L) %in% full$selected))       # the exact copies are pruned

  # the early-stopping walk yields precisely the K most relevant survivors
  for (k in seq_len(length(full$selected))) {
    expect_identical(fcbf_select(X, y, delta = 0, top = k)$selected,
                     utils::head(full$selected, k))
  }
  # `top` is a cap: never exceeded, and capping above the survivor count is harmless
  expect_lte(length(fcbf_select(X, y, delta = 0, top = 1)$selected), 1L)
  expect_identical(fcbf_select(X, y, delta = 0, top = 99L)$selected, full$selected)

  # the trace stays well-formed under a cap (only examined features are traced)
  tk <- fcbf_select(X, y, delta = 0, top = 2)
  expect_true(all(tk$trace$removed$su_blanket_feature >= tk$trace$removed$su_feature_target))
  expect_named(tk$trace$representative_per_cluster, as.character(tk$selected))
})

test_that("fcbf_select rejects an invalid top", {
  X <- matrix(c(0, 1, 0, 1, 1, 0), nrow = 3)
  expect_error(fcbf_select(X, y = c(0, 1, 1), top = 0), "positive integer")
  expect_error(fcbf_select(X, y = c(0, 1, 1), top = -2), "positive integer")
  expect_error(fcbf_select(X, y = c(0, 1, 1), top = c(1, 2)), "positive integer")
})
