# Acceptance tests for make_sigma_blocks().
#
# A small but structurally complete configuration (mirrors the ficha shape at
# reduced size so eigen() stays cheap):
#   clusters:    k_drivers * (cluster_size + 1) = 4 * (3 + 1) = 16 columns
#   distractors: n_distractor_subblocks * distractor_size = 3 * 4 = 12 columns
#   interaction: 2 columns
#   noise:       p - 30 = 60 - 30 = 30 columns
make_cfg <- function(p = 60, k_drivers = 4, cluster_size = 3, rho = 0.8,
                     n_distractor_subblocks = 3, distractor_size = 4,
                     distractor_rho = 0.7, interaction = TRUE) {
  list(p = p, k_drivers = k_drivers, cluster_size = cluster_size, rho = rho,
       n_distractor_subblocks = n_distractor_subblocks,
       distractor_size = distractor_size, distractor_rho = distractor_rho,
       interaction = interaction)
}

build <- function(cfg) do.call(make_sigma_blocks, cfg)

test_that("Sigma is square p x p, symmetric, unit-diagonal and positive definite", {
  cfg <- make_cfg()
  Sigma <- build(cfg)$Sigma

  expect_true(is.matrix(Sigma))
  expect_equal(dim(Sigma), c(cfg$p, cfg$p))
  expect_true(isSymmetric(Sigma))
  expect_equal(diag(Sigma), rep(1, cfg$p))

  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("manifest role counts match the arguments", {
  cfg <- make_cfg()
  m <- build(cfg)$manifest

  expect_equal(nrow(m), cfg$p)
  expect_equal(sum(m$role == "driver"),      cfg$k_drivers)
  expect_equal(sum(m$role == "redundant"),   cfg$k_drivers * cfg$cluster_size)
  expect_equal(sum(m$role == "distractor"),  cfg$n_distractor_subblocks * cfg$distractor_size)
  expect_equal(sum(m$role == "interaction"), 2L)
  expected_noise <- cfg$p -
    cfg$k_drivers * (cfg$cluster_size + 1) -
    cfg$n_distractor_subblocks * cfg$distractor_size - 2L
  expect_equal(sum(m$role == "noise"), expected_noise)
})

test_that("cluster_id groups driver clusters and distractor subblocks correctly", {
  cfg <- make_cfg()
  m <- build(cfg)$manifest

  # one cluster id per driver, each holding exactly 1 driver + cluster_size redundants
  driver_ids <- m$cluster_id[m$role == "driver"]
  expect_equal(length(unique(driver_ids)), cfg$k_drivers)
  for (cid in driver_ids) {
    members <- m[m$cluster_id == cid & m$role %in% c("driver", "redundant"), ]
    expect_equal(nrow(members), cfg$cluster_size + 1L)
    expect_equal(sum(members$role == "driver"), 1L)
    expect_equal(sum(members$role == "redundant"), cfg$cluster_size)
  }

  # one cluster id per distractor subblock, each of size distractor_size
  distractor_ids <- m$cluster_id[m$role == "distractor"]
  expect_equal(length(unique(distractor_ids)), cfg$n_distractor_subblocks)
  expect_true(all(table(distractor_ids) == cfg$distractor_size))

  # driver and distractor id spaces do not overlap
  expect_length(intersect(unique(driver_ids), unique(distractor_ids)), 0)

  # noise has no cluster membership
  expect_true(all(is.na(m$cluster_id[m$role == "noise"])))
})

test_that("driver-driver latent correlation is zero", {
  cfg <- make_cfg()
  res <- build(cfg)
  driver_cols <- which(res$manifest$role == "driver")
  sub <- res$Sigma[driver_cols, driver_cols]
  expect_equal(sub, diag(length(driver_cols)))
})

test_that("intra-cluster correlation equals rho and distractor blocks equal distractor_rho", {
  cfg <- make_cfg()
  res <- build(cfg)
  m <- res$manifest

  for (cid in unique(m$cluster_id[m$role == "driver"])) {
    idx <- which(m$cluster_id == cid & m$role %in% c("driver", "redundant"))
    block <- res$Sigma[idx, idx]
    expect_true(all(abs(block[upper.tri(block)] - cfg$rho) < 1e-12))
  }

  for (cid in unique(m$cluster_id[m$role == "distractor"])) {
    idx <- which(m$cluster_id == cid)
    block <- res$Sigma[idx, idx]
    expect_true(all(abs(block[upper.tri(block)] - cfg$distractor_rho) < 1e-12))
  }
})

test_that("only drivers carry beta and driver betas are mixed sign", {
  cfg <- make_cfg()
  m <- build(cfg)$manifest

  driver_beta <- m$beta[m$role == "driver"]
  expect_true(all(driver_beta != 0))
  expect_true(any(driver_beta > 0) && any(driver_beta < 0))   # mixed signs
  expect_true(all(m$beta[m$role != "driver"] == 0))
})

test_that("interaction = FALSE drops the interaction pair and stays positive definite", {
  cfg <- make_cfg(interaction = FALSE)
  res <- build(cfg)
  m <- res$manifest

  expect_equal(sum(m$role == "interaction"), 0L)
  expected_noise <- cfg$p -
    cfg$k_drivers * (cfg$cluster_size + 1) -
    cfg$n_distractor_subblocks * cfg$distractor_size
  expect_equal(sum(m$role == "noise"), expected_noise)

  ev <- eigen(res$Sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("guard: chosen rho keeps every compound-symmetry block positive definite", {
  cfg <- make_cfg()
  res <- build(cfg)
  m <- res$manifest

  group_ids <- unique(m$cluster_id[m$role %in% c("driver", "redundant", "distractor")])
  for (cid in group_ids) {
    idx <- which(m$cluster_id == cid)
    block <- res$Sigma[idx, idx]
    ev <- eigen(block, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(ev > 0))
  }
})

test_that("guard: an out-of-range rho or distractor_rho is rejected", {
  # rho = 1 makes the cluster block singular (1 - rho = 0)
  expect_error(
    make_sigma_blocks(p = 60, k_drivers = 4, cluster_size = 3, rho = 1,
                      n_distractor_subblocks = 3, distractor_size = 4,
                      distractor_rho = 0.7),
    "positive definite"
  )
  # distractor_rho = 1 likewise
  expect_error(
    make_sigma_blocks(p = 60, k_drivers = 4, cluster_size = 3, rho = 0.8,
                      n_distractor_subblocks = 3, distractor_size = 4,
                      distractor_rho = 1),
    "positive definite"
  )
  # cluster block size = cluster_size + 1 = 4, so rho must exceed -1/3
  expect_error(
    make_sigma_blocks(p = 60, k_drivers = 4, cluster_size = 3, rho = -0.5,
                      n_distractor_subblocks = 3, distractor_size = 4,
                      distractor_rho = 0.7),
    "positive definite"
  )
})

test_that("p too small to hold the structured columns is an error", {
  expect_error(
    make_sigma_blocks(p = 10, k_drivers = 4, cluster_size = 3, rho = 0.8,
                      n_distractor_subblocks = 3, distractor_size = 4,
                      distractor_rho = 0.7),
    "too small"
  )
})

test_that("the default ficha dimensions hold (8 drivers + 40 redundants = 48 columns)", {
  res <- make_sigma_blocks(p = 5000, k_drivers = 8, cluster_size = 5, rho = 0.8,
                           n_distractor_subblocks = 20, distractor_size = 10,
                           distractor_rho = 0.7, interaction = TRUE)
  m <- res$manifest
  expect_equal(nrow(m), 5000L)
  expect_equal(sum(m$role == "driver"), 8L)
  expect_equal(sum(m$role == "redundant"), 40L)
  expect_equal(sum(m$role %in% c("driver", "redundant")), 48L)
  expect_equal(sum(m$role == "distractor"), 200L)
  expect_equal(sum(m$role == "interaction"), 2L)
  expect_equal(sum(m$role == "noise"), 4750L)
})
