# Validation of fcbf_select() against the reference implementation, the FCBF package
# (Bioconductor until release 3.19). The expected selections below are the output of
# FCBF::fcbf() version 1.99.1 (github.com/lubianat/FCBF, commit 51b9600, the version pinned
# in renv.lock), computed once and frozen here, so the check runs without that package. The
# full comparison, on 400 generated datasets and a single-cell dataset, is
# experiments/benchmark/realdata/fcbf_crossvalidation.R.

test_that("fcbf_select selects the same set as the reference FCBF on non-degenerate data", {
  set.seed(11)
  n <- 120
  y <- rep(1:3, each = 40)
  noisy <- function(v, k) { v[sample(length(v), k)] <- sample(1:3, k, TRUE); v }
  A <- noisy(y, 12); Ared <- noisy(A, 6); B <- noisy(y, 24)
  N1 <- rep_len(c(1, 2, 3, 1), n); N2 <- rep_len(c(2, 1), n)
  X <- cbind(A = A, Ared = Ared, B = B, N1 = N1, N2 = N2)

  ## FCBF::fcbf(minimum_su = delta) keeps A and B at every threshold: Ared is pruned as
  ## redundant with A, and N1 and N2 carry no signal
  reference <- c(1L, 3L)
  for (d in c(0.05, 0.15, 0.30)) {
    expect_setequal(fcbf_select(X, y, delta = d)$selected, reference)
  }
})

test_that("the >= convention prunes a target-identical duplicate that FCBF's > keeps", {
  y <- rep(c(1, 2, 3), each = 15)
  f <- y                                  # feature identical to the target: SU(f, y) = 1
  X <- cbind(f1 = f, f2 = f)              # exact duplicate, both equal the target

  ## FCBF::fcbf(minimum_su = 0.05) keeps both columns: its strict > does not prune f2 when
  ## SU(f1, f2) = SU(f2, y) = 1, the only boundary where the two implementations differ
  reference <- c(1L, 2L)
  mine <- fcbf_select(X, y, delta = 0.05)$selected
  expect_equal(mine, 1L)                  # the published >= rule prunes the duplicate
  expect_equal(setdiff(reference, mine), 2L)  # and that duplicate is the only difference
})
