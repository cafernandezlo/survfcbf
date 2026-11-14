# Acceptance tests for score_selection_clusterwise().
# Layout for p = 30, k_drivers = 3, cluster_size = 3, 2 distractor subblocks x 4:
#   clusters cid 1/2/3 at cols 1-4 / 5-8 / 9-12 (driver = first col of each block),
#   distractors cid 4/5 at cols 13-16 / 17-20, interaction cols 21-22, noise 23-30.
man <- function() {
  make_sigma_blocks(p = 30, k_drivers = 3, cluster_size = 3, rho = 0.8,
                    n_distractor_subblocks = 2, distractor_size = 4,
                    distractor_rho = 0.7, interaction = TRUE)$manifest
}

test_that("selecting one driver per cluster is perfect, parsimonious recovery", {
  m <- man()
  drivers <- which(m$role == "driver")          # 1, 5, 9
  s <- score_selection_clusterwise(drivers, m)

  expect_equal(s$sensitivity, 1)
  expect_equal(s$fdr, 0)
  expect_equal(s$f1, 1)
  expect_true(s$exact_recovery)
  expect_equal(s$parsimony_penalty, 0)
  expect_equal(unname(s$representative_per_cluster[["1"]]), drivers[1])
  expect_length(s$representative_per_cluster, 3L)
})

test_that("two members of one cluster raise parsimony but not sensitivity", {
  m <- man()
  drivers <- which(m$role == "driver")
  red1 <- which(m$role == "redundant" & m$cluster_id == m$cluster_id[drivers[1]])[1]
  s <- score_selection_clusterwise(c(drivers, red1), m)

  expect_equal(s$sensitivity, 1)               # still all clusters recovered
  expect_equal(s$parsimony_penalty, 1)         # one redundant extra in cluster 1
  expect_true(s$exact_recovery)                # parsimony does not break exact recovery
  # driver listed first in the cluster's representatives
  expect_equal(s$representative_per_cluster[[as.character(m$cluster_id[drivers[1]])]][1], drivers[1])
})

test_that("a distractor selection is a false discovery and breaks exact recovery", {
  m <- man()
  d1 <- which(m$role == "driver")[1]
  distractor <- which(m$role == "distractor")[1]
  s <- score_selection_clusterwise(c(d1, distractor), m)

  expect_equal(s$sensitivity, 1 / 3)
  expect_equal(s$fdr, 0.5)                      # 1 of 2 selections is a distractor
  expect_equal(s$f1, 0.4)                       # P = 1/2, R = 1/3
  expect_false(s$exact_recovery)
  expect_equal(s$parsimony_penalty, 0)
})

test_that("interaction selections are neither true positives nor false discoveries", {
  m <- man()
  drivers <- which(m$role == "driver")
  inter <- which(m$role == "interaction")[1]
  s <- score_selection_clusterwise(c(drivers, inter), m)

  expect_equal(s$sensitivity, 1)
  expect_equal(s$fdr, 0)                        # interaction not counted as a false discovery
  expect_true(s$exact_recovery)                # and does not break exact recovery
  expect_equal(s$f1, 1)
})

test_that("the empty selection scores zero everywhere", {
  m <- man()
  s <- score_selection_clusterwise(integer(0), m)
  expect_equal(s$sensitivity, 0)
  expect_equal(s$fdr, 0)
  expect_equal(s$f1, 0)
  expect_false(s$exact_recovery)
  expect_equal(s$parsimony_penalty, 0)
  expect_length(s$representative_per_cluster, 0L)
})

test_that("logical and integer selections agree, and bad indices are rejected", {
  m <- man()
  drivers <- which(m$role == "driver")
  lg <- logical(nrow(m)); lg[drivers] <- TRUE
  expect_equal(score_selection_clusterwise(lg, m), score_selection_clusterwise(drivers, m))

  expect_error(score_selection_clusterwise(c(1L, 999L), m), "1..nrow")
  expect_error(score_selection_clusterwise(c(TRUE, FALSE), m), "length nrow")
})
