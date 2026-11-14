# Acceptance tests for survfcbf_select().

# Small, fast scenario for the structural tests (reduced distractor block).
small_sim <- function(seed) {
  simulate_dataset(default_params(
    p = 120, k_drivers = 3, cluster_size = 3,
    n_distractor_subblocks = 3, distractor_size = 4,
    n = 300, seed = seed, calib_n_presample = 1000))
}

# Silences only the follow-up warning, in tests that run heavily censored data for
# another purpose. Any other warning still reaches testthat.
muffle_followup <- function(expr) {
  withCallingHandlers(expr, survfcbf_followup_warning = function(w) invokeRestart("muffleWarning"))
}

test_that("survfcbf_select recovers concentrated driver signal with low FDR", {
  # concentrated signal (few strong drivers, high C, light censoring) -> clean recovery
  sim <- simulate_dataset(default_params(
    p = 300, k_drivers = 3, cluster_size = 5, n = 800,
    censor_fraction = 0.10, target_C = 0.90, seed = 31, calib_n_presample = 1500))
  res <- survfcbf_select(sim$time, sim$status, sim$X, q = 10, B = 50, seed = 32)
  sc <- score_selection_clusterwise(res$selected, sim$manifest)

  expect_gte(sc$sensitivity, 2 / 3)     # recovers at least 2 of the 3 driver clusters
  expect_lte(sc$fdr, 0.2)               # essentially no distractor/noise selections
  expect_gte(length(res$selected), 2L)
})

test_that("the output carries everything needed to apply to held-out data", {
  sim <- small_sim(3)
  res <- survfcbf_select(sim$time, sim$status, sim$X, q = 8, B = 30, seed = 4)

  expect_setequal(names(res), c("selected", "cutpoints", "delta", "target", "relevance", "trace"))
  expect_length(res$cutpoints, ncol(sim$X))
  expect_length(res$target, nrow(sim$X))
  expect_length(res$relevance, ncol(sim$X))
  expect_true(res$delta >= 0 && res$delta <= 1)

  # training cutpoints re-bin fresh data with no leakage
  j <- if (length(res$selected)) res$selected[1] else 1L
  rb <- discretize_quantile(stats::rnorm(40), breaks = res$cutpoints[[j]])$binned
  expect_false(anyNA(rb))
})

test_that("a fixed delta makes the selection deterministic", {
  sim <- small_sim(5)
  a <- survfcbf_select(sim$time, sim$status, sim$X, q = 8, delta = 0.1)
  b <- survfcbf_select(sim$time, sim$status, sim$X, q = 8, delta = 0.1)
  expect_identical(a$selected, b$selected)
  expect_equal(a$delta, 0.1)
})

test_that("survfcbf_select is reproducible under a fixed seed", {
  sim <- small_sim(6)
  a <- survfcbf_select(sim$time, sim$status, sim$X, q = 8, B = 30, seed = 7)
  b <- survfcbf_select(sim$time, sim$status, sim$X, q = 8, B = 30, seed = 7)
  expect_identical(a$selected, b$selected)
  expect_equal(a$delta, b$delta)
})

test_that("cox_score relevance recovers concentrated signal and stays SU-comparable", {
  skip_if_not_installed("survival")
  sim <- simulate_dataset(default_params(
    p = 300, k_drivers = 3, cluster_size = 5, n = 800,
    censor_fraction = 0.10, target_C = 0.90, seed = 31, calib_n_presample = 1500))
  res <- survfcbf_select(sim$time, sim$status, sim$X, q = 10,
                         relevance = "cox_score", delta = 0.15)

  expect_true(all(res$relevance >= 0 & res$relevance <= 1))   # rescaled to [0, 1]
  expect_null(res$target)                                     # no martingale target in this ablation
  sc <- score_selection_clusterwise(res$selected, sim$manifest)
  expect_gte(sc$sensitivity, 2 / 3)
  expect_lte(sc$fdr, 0.3)
})

test_that("cox_score relevance requires an explicit delta", {
  sim <- small_sim(8)
  expect_error(
    survfcbf_select(sim$time, sim$status, sim$X, relevance = "cox_score"),
    "explicit"
  )
})

test_that("a top cap selects an explicit cardinality, deterministically and without a permutation null", {
  sim <- small_sim(11)
  a <- survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 4)   # no seed needed
  b <- survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 4)
  expect_identical(a$selected, b$selected)                            # deterministic: no random delta null
  expect_gte(length(a$selected), 2L)
  expect_lte(length(a$selected), 4L)                                  # cardinality capped at top
  expect_setequal(names(a), c("selected", "cutpoints", "delta", "target", "relevance", "trace"))
  expect_null(a$delta)                                                # no threshold used; cardinality by top
  expect_length(a$cutpoints, ncol(sim$X))
  expect_length(a$target, nrow(sim$X))
})

test_that("top selection stays non-empty and parsimonious in a low-event, high-noise regime", {
  # the genomic regime that collapses the permutation-null delta self-selection
  sim <- simulate_dataset(default_params(
    p = 300, k_drivers = 3, cluster_size = 5, n = 300,
    censor_fraction = 0.75, target_C = 0.70, seed = 21, calib_n_presample = 1000))
  res <- muffle_followup(survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 5))
  expect_gte(length(res$selected), 1L)                                # selects a set where delta mode would not
  expect_lte(length(res$selected), 5L)
})

test_that("a top cap lifts the cox_score explicit-delta requirement", {
  skip_if_not_installed("survival")
  sim <- small_sim(12)
  res <- survfcbf_select(sim$time, sim$status, sim$X, relevance = "cox_score", q = 2, top = 4)
  expect_lte(length(res$selected), 4L)
  expect_null(res$target)                                             # cox_score keeps no martingale target
  expect_true(all(res$relevance >= 0 & res$relevance <= 1))
})

test_that("survfcbf_select rejects an invalid top", {
  sim <- small_sim(13)
  expect_error(survfcbf_select(sim$time, sim$status, sim$X, top = 0), "positive integer")
  expect_error(survfcbf_select(sim$time, sim$status, sim$X, top = -1), "positive integer")
})

test_that("ipcw_weights uses the left-continuous G(t-) at tied times", {
  skip_if_not_installed("survival")
  # one censoring at t = 2; the event tying it must use G(2-) = 1, not G(2) = 2/3
  time   <- c(1, 2, 2, 3)
  status <- c(1, 0, 1, 1)
  w <- ipcw_weights(time, status)
  expect_equal(w, c(1, 0, 1, 1.5))   # G(t-) = 1, 1, 1, 2/3 -> weights 1, 0, 1, 1.5
})

test_that("ipcw_weights rejects all-censored input", {
  expect_error(ipcw_weights(c(1, 2, 3), c(0, 0, 0)), "at least one event")
})

test_that("the ipcw_binned target runs end to end with IPCW-weighted SU", {
  skip_if_not_installed("survival")
  sim <- simulate_dataset(default_params(
    p = 300, k_drivers = 3, cluster_size = 5, n = 800,
    censor_fraction = 0.30, target_C = 0.90, seed = 31, calib_n_presample = 1500))
  res <- survfcbf_select(sim$time, sim$status, sim$X, q = 10, B = 50,
                         target_method = "ipcw_binned", seed = 32)

  expect_true(all(res$relevance >= 0 & res$relevance <= 1))   # weighted SU stays in [0, 1]
  expect_equal(res$target, sim$time)                          # ipcw_binned target is the observed time
  sc <- score_selection_clusterwise(res$selected, sim$manifest)
  expect_gte(sc$sensitivity, 1 / 3)                           # recovers signal (weaker leg under censoring)
  expect_lte(sc$fdr, 0.4)
})

## Column-order invariance. At q = 2 SU takes few distinct values, so many features tie
## on relevance and, without a secondary key, the column order decided which of them
## entered the panel.
order_sim <- function() simulate_dataset(default_params(
  p = 400, k_drivers = 3, cluster_size = 5, n_distractor_subblocks = 4, distractor_size = 8,
  n = 200, censor_fraction = 0.7, seed = 41, calib_n_presample = 800))

test_that("relevance at q = 2 is heavily tied, which is why a tie-break is needed", {
  sim <- order_sim()
  res <- survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 10)
  expect_lt(length(unique(round(res$relevance, 12))), ncol(sim$X) / 10)
})

test_that("the selection does not depend on the column order of X", {
  sim <- order_sim()
  set.seed(42)
  perm <- sample(ncol(sim$X))
  for (args in list(list(q = 2, top = 10), list(q = 2, delta = 0.005), list(q = 4, top = 10))) {
    a <- do.call(survfcbf_select, c(list(sim$time, sim$status, sim$X), args))
    b <- do.call(survfcbf_select, c(list(sim$time, sim$status, sim$X[, perm]), args))
    expect_identical(perm[b$selected], a$selected)            # same genes, in the same order
  }
})

test_that("tiebreak = 'index' reproduces the column-order rule of fcbf_select", {
  sim <- order_sim()
  res <- survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 8, tiebreak = "index")
  ym <- survival_target(sim$time, sim$status)
  Xd <- vapply(seq_len(ncol(sim$X)), function(j) discretize_quantile(sim$X[, j], q = 2)$binned,
               integer(nrow(sim$X)))
  yd <- discretize_quantile(ym, q = 2)$binned
  expect_identical(res$selected, fcbf_select(Xd, yd, top = 8)$selected)
})

test_that("the Spearman tie-break is computed on the supplied data only", {
  sim <- order_sim()
  tb <- spearman_tiebreak(sim$X, survival_target(sim$time, sim$status))
  expect_length(tb, ncol(sim$X))
  expect_true(all(tb >= 0 & tb <= 1))
  expect_identical(spearman_tiebreak(cbind(sim$X[, 1:3], 1), sim$time)[4], 0)   # constant feature
})

## Follow-up warning. With heavy censoring the top bin of the martingale target holds the
## events plus the censored patients with the shortest follow-up, so relevance can track
## follow-up rather than risk. The warning reports it and must never change the selection.
followup_sim <- function(censor_fraction) simulate_dataset(default_params(
  p = 120, k_drivers = 3, cluster_size = 3, n_distractor_subblocks = 3, distractor_size = 4,
  n = 300, censor_fraction = censor_fraction, seed = 51, calib_n_presample = 1000))

test_that("a censored majority in the top target bin raises a classed follow-up warning", {
  sim <- followup_sim(0.85)
  tr <- 1:220                                                  # a training fold: the share is computed on it
  yd <- discretize_quantile(survival_target(sim$time[tr], sim$status[tr]), q = 2)$binned
  top <- yd == max(yd)
  share <- mean(sim$status[tr][top] == 0)
  expect_gt(share, 0.5)

  w <- tryCatch(survfcbf_select(sim$time[tr], sim$status[tr], sim$X[tr, ], q = 2, top = 5),
                survfcbf_followup_warning = function(w) w)
  expect_s3_class(w, c("survfcbf_followup_warning", "warning", "condition"), exact = TRUE)
  expect_equal(w$prop, share)
  expect_identical(w$n_top, sum(top))
  expect_match(conditionMessage(w), sprintf("%.1f%% of the patients in the top bin", 100 * share), fixed = TRUE)
  expect_match(conditionMessage(w), "may partly reflect follow-up length rather than risk", fixed = TRUE)
  expect_match(conditionMessage(w), "larger q", fixed = TRUE)
})

test_that("the follow-up warning fires in top-K and threshold mode and for any q", {
  sim <- followup_sim(0.92)
  for (args in list(list(q = 2, top = 5), list(q = 2, delta = 0.01), list(q = 2, B = 20, seed = 1),
                    list(q = 4, top = 5), list(q = 4, delta = 0.01))) {
    expect_warning(do.call(survfcbf_select, c(list(sim$time, sim$status, sim$X), args)),
                   class = "survfcbf_followup_warning")
  }
})

test_that("light censoring, or a q whose top bin holds only events, raises no warning", {
  sim <- followup_sim(0.10)
  expect_no_warning(survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 5))
  expect_no_warning(survfcbf_select(sim$time, sim$status, sim$X, q = 2, delta = 0.01))

  # with 85% censored the top decile is all events: the check follows the chosen q
  sim <- followup_sim(0.85)
  yd <- discretize_quantile(survival_target(sim$time, sim$status), q = 10)$binned
  expect_lte(mean(sim$status[yd == max(yd)] == 0), 0.5)
  expect_no_warning(survfcbf_select(sim$time, sim$status, sim$X, q = 10, top = 5),
                    class = "survfcbf_followup_warning")
})

test_that("the follow-up warning leaves the selection untouched", {
  sim <- followup_sim(0.85)
  for (args in list(list(q = 2, top = 5), list(q = 2, delta = 0.01), list(q = 2, B = 20, seed = 3),
                    list(q = 2, B = 20), list(q = 2, top = 5, tiebreak = "index"))) {
    run <- function() do.call(survfcbf_select, c(list(sim$time, sim$status, sim$X), args))
    set.seed(9)                                                # seed = NULL draws from the global stream
    expect_warning(loud <- run(), class = "survfcbf_followup_warning")
    set.seed(9)
    quiet <- muffle_followup(run())
    expect_identical(quiet, loud)
  }
})

test_that("the follow-up check covers only the martingale target with SU relevance", {
  sim <- followup_sim(0.92)
  expect_no_warning(survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 5, relevance = "cox_score"),
                    class = "survfcbf_followup_warning")
  expect_no_warning(survfcbf_select(sim$time, sim$status, sim$X, q = 2, delta = 0.01,
                                    target_method = "ipcw_binned"),
                    class = "survfcbf_followup_warning")
})
