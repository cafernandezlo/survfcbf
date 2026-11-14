# Acceptance tests for the baseline interface and sanity-floor selectors.

test_that("cut_ranking returns the top features by score, ties by index", {
  sc <- c(0.1, 0.9, 0.5, 0.9)
  expect_equal(cut_ranking(sc, k = 1), 2L)              # single top
  expect_equal(cut_ranking(sc, k = 2), c(2L, 4L))       # tie at 0.9 -> lower index first
  expect_equal(sort(cut_ranking(sc, prop = 0.5)), c(2L, 4L))   # round(0.5*4) = 2
  expect_length(cut_ranking(sc, prop = 1), 4L)          # keep all
  expect_equal(cut_ranking(sc, k = 99), c(2L, 4L, 3L, 1L))     # k capped at p, full order
  expect_error(cut_ranking(sc), "either")               # needs prop or k
  expect_error(cut_ranking(c(1, NA, 2), k = 1), "NA")
})

test_that("variance_select ranks by marginal variance and cuts at prop", {
  set.seed(1)
  X <- cbind(stats::rnorm(200, sd = 5), stats::rnorm(200, sd = 1), stats::rnorm(200, sd = 3))
  tt <- stats::rexp(200); st <- rep(1L, 200)
  r <- variance_select(tt, st, X, prop = 1 / 3)
  expect_equal(which.max(r$score), 1L)                  # highest-variance column scores highest
  expect_equal(r$selected, 1L)                          # top 1 of 3
  expect_null(variance_select(tt, st, X)$selected)      # no prop -> pure ranker
  expect_length(r$score, 3L)
})

test_that("nofilter_select keeps every feature", {
  X <- matrix(stats::rnorm(50), 10, 5)
  r <- nofilter_select(rep(1, 10), rep(1L, 10), X)
  expect_equal(r$selected, 1:5)
  expect_length(r$score, 5L)
})

test_that("sanity-floor selectors integrate with cluster-level scoring on simulated data", {
  sim <- simulate_dataset(default_params(
    p = 200, k_drivers = 3, cluster_size = 4,
    n_distractor_subblocks = 3, distractor_size = 4,
    n = 400, seed = 5, calib_n_presample = 1000))

  vs <- variance_select(sim$time, sim$status, sim$X, prop = 0.1)
  expect_length(vs$selected, max(1L, round(0.1 * 200)))
  sc <- score_selection_clusterwise(vs$selected, sim$manifest)
  expect_true(all(c("sensitivity", "fdr") %in% names(sc)))
  expect_true(sc$fdr >= 0 && sc$fdr <= 1)

  nf <- nofilter_select(sim$time, sim$status, sim$X)
  scn <- score_selection_clusterwise(nf$selected, sim$manifest)
  expect_equal(scn$sensitivity, 1)                      # keeping all features recovers every cluster
})

test_that("baseline_registry lists the wrapped selectors with valid metadata", {
  reg <- baseline_registry()
  expect_gte(length(reg), 9L)
  nms <- vapply(reg, `[[`, character(1), "name")
  expect_setequal(nms, c("variance", "nofilter", "mim_surv", "mrmr_surv", "carss",
                         "cox_glmnet_l1", "cox_glmnet_l2", "rsf", "grplasso_ensemble"))
  expect_true(all(vapply(reg, function(e) is.function(e$fn), logical(1))))
  expect_true(all(vapply(reg, function(e) e$type %in% c("set", "ranker"), logical(1))))
  expect_true(all(vapply(reg, function(e) is.logical(e$self_selects), logical(1))))
  carss_e <- Filter(function(e) e$name == "carss", reg)[[1]]   # ranked at prop (Bommert), not self-selected
  expect_identical(carss_e$type, "ranker")
  expect_false(carss_e$self_selects)
})

test_that("run_baseline invokes an entry and threads prop, seed and fixed args", {
  sim <- simulate_dataset(default_params(
    p = 150, k_drivers = 3, cluster_size = 4, n = 400,
    n_distractor_subblocks = 4, distractor_size = 6,
    seed = 9, calib_n_presample = 1000))
  reg <- baseline_registry()
  byname <- function(n) Filter(function(e) e$name == n, reg)[[1]]

  # base-only entries run with no Suggested package present
  v <- run_baseline(byname("variance"), sim$time, sim$status, sim$X, prop = 0.1)
  expect_length(v$score, 150L)
  expect_length(v$selected, max(1L, round(0.1 * 150)))
  expect_length(run_baseline(byname("nofilter"), sim$time, sim$status, sim$X)$selected, 150L)

  # cox_glmnet_l1 carries alpha = 1 through args -> the nonzero set
  skip_if_not_installed("glmnet")
  l1 <- run_baseline(byname("cox_glmnet_l1"), sim$time, sim$status, sim$X, seed = 1)
  expect_setequal(l1$selected, which(l1$coef != 0))
})
