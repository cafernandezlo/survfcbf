# Acceptance tests for sample_features().

sigma_small <- function() {
  make_sigma_blocks(p = 24, k_drivers = 3, cluster_size = 3, rho = 0.8,
                    n_distractor_subblocks = 2, distractor_size = 4,
                    distractor_rho = 0.7, interaction = TRUE)
}

test_that("gaussian marginal reproduces Sigma within Monte Carlo error", {
  bl <- sigma_small()
  res <- sample_features(n = 5000, Sigma = bl$Sigma, marginals = "gaussian", seed = 42)

  expect_equal(dim(res$X), c(5000L, 24L))
  expect_false(is.null(res$realized_cor))
  expect_true(max(abs(res$realized_cor - bl$Sigma)) < 0.08)

  m <- bl$manifest
  cl <- which(m$cluster_id == 1 & m$role %in% c("driver", "redundant"))
  block <- res$realized_cor[cl, cl]
  expect_true(all(abs(block[upper.tri(block)] - 0.8) < 0.06))   # intra-cluster ~ rho

  dcols <- which(m$role == "driver")
  dd <- res$realized_cor[dcols, dcols]
  expect_true(all(abs(dd[upper.tri(dd)]) < 0.06))               # driver-driver ~ 0
})

test_that("negbin marginal matches requested moments and is count-valued", {
  bl <- sigma_small()
  res <- sample_features(n = 6000, Sigma = bl$Sigma, marginals = "negbin",
                         negbin_mu = 50, negbin_size = 2, seed = 7, return_cor = FALSE)
  X <- res$X

  expect_null(res$realized_cor)
  expect_true(all(X >= 0))
  expect_true(all(X == round(X)))                               # integer counts

  expect_true(abs(mean(colMeans(X)) - 50) < 2)                  # mean ~ mu
  theo_var <- 50 + 50^2 / 2                                     # negbin var = mu + mu^2/size
  expect_true(abs(mean(apply(X, 2L, stats::var)) / theo_var - 1) < 0.15)
})

test_that("copula attenuates the realised intra-cluster correlation of negbin margins", {
  bl <- sigma_small()
  res <- sample_features(n = 6000, Sigma = bl$Sigma, marginals = "negbin",
                         negbin_mu = 50, negbin_size = 2, seed = 11)
  m <- bl$manifest
  cl <- which(m$cluster_id == 1 & m$role %in% c("driver", "redundant"))
  off <- res$realized_cor[cl, cl][upper.tri(diag(length(cl)))]
  expect_true(all(off > 0.2))                                   # still substantially correlated
  expect_true(all(off < 0.8))                                   # but attenuated below latent rho
})

test_that("empirical marginal reproduces the assigned reference gene distribution", {
  n_genes <- 40L; n_samp <- 200L
  spread <- stats::qnorm(stats::ppoints(n_samp))                # fixed, mean 0
  ref <- t(sapply(seq_len(n_genes), function(g) 100 + 5 * g + spread))
  rownames(ref) <- paste0("g", seq_len(n_genes))

  bl <- sigma_small()                                           # p = 24
  er <- build_empirical_reference(ref, p = ncol(bl$Sigma), seed = 3)
  res <- sample_features(n = 6000, Sigma = bl$Sigma, marginals = "empirical",
                         reference = er, seed = 5, return_cor = FALSE)

  assigned_g <- as.integer(sub("^g", "", er$genes))
  expected_mean <- 100 + 5 * assigned_g
  expect_equal(colMeans(res$X), expected_mean, tolerance = 0.02)
})

test_that("sampling is reproducible and honours guards", {
  bl <- sigma_small()
  a <- sample_features(n = 200, Sigma = bl$Sigma, marginals = "gaussian", seed = 99)
  b <- sample_features(n = 200, Sigma = bl$Sigma, marginals = "gaussian", seed = 99)
  expect_identical(a$X, b$X)

  expect_null(sample_features(n = 50, Sigma = bl$Sigma, marginals = "gaussian",
                              return_cor = FALSE)$realized_cor)

  expect_error(
    sample_features(n = 50, Sigma = bl$Sigma, marginals = "empirical"),
    "reference"
  )
  ref <- matrix(seq_len(200), nrow = 10)                        # 10 genes x 20 samples
  er <- build_empirical_reference(ref, p = 5, seed = 1)         # p = 5 != 24
  expect_error(
    sample_features(n = 50, Sigma = bl$Sigma, marginals = "empirical", reference = er),
    "built for p"
  )
})

test_that("block-aware latent sampling recovers a known block covariance and handles n = 1", {
  bl <- make_sigma_blocks(p = 20, k_drivers = 2, cluster_size = 3, rho = 0.7,
                          n_distractor_subblocks = 1, distractor_size = 4,
                          distractor_rho = 0.5, interaction = TRUE)
  res <- sample_features(n = 20000, Sigma = bl$Sigma, marginals = "gaussian", seed = 123)
  expect_lt(max(abs(res$realized_cor - bl$Sigma)), 0.06)      # every block recovered, zeros between

  one <- sample_features(n = 1, Sigma = bl$Sigma, marginals = "gaussian", seed = 1, return_cor = FALSE)
  expect_equal(dim(one$X), c(1L, 20L))                        # n = 1 stays a 1 x p matrix
})

