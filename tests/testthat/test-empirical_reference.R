# Acceptance tests for build_empirical_reference().

# Deterministic synthetic reference: gene g has mean 10 + g and a fixed spread,
# all values positive so every gene is "expressed".
make_ref <- function(n_genes = 40L, n_samp = 200L) {
  spread <- stats::qnorm(stats::ppoints(n_samp))      # mean 0, sd ~ 1, fixed
  ref <- t(sapply(seq_len(n_genes), function(g) 10 + g + spread))
  rownames(ref) <- paste0("g", seq_len(n_genes))
  ref
}

test_that("reference has the documented structure and sorted quantile columns", {
  ref <- make_ref()
  er <- build_empirical_reference(ref, p = 12, seed = 1)

  expect_s3_class(er, "survfcbf_empirical_reference")
  expect_equal(er$p, 12L)
  expect_equal(er$n_samples, 200L)
  expect_equal(er$n_pool, 40L)
  expect_equal(dim(er$quantiles), c(200L, 12L))
  expect_length(er$genes, 12L)
  expect_true(all(er$genes %in% rownames(ref)))

  for (j in seq_len(ncol(er$quantiles))) {
    expect_false(is.unsorted(er$quantiles[, j]))       # each column sorted ascending
  }
})

test_that("gene-to-feature assignment is reproducible under a fixed seed", {
  ref <- make_ref()
  a <- build_empirical_reference(ref, p = 15, seed = 123)
  b <- build_empirical_reference(ref, p = 15, seed = 123)
  expect_identical(a$genes, b$genes)
  expect_identical(a$quantiles, b$quantiles)
})

test_that("degenerate genes (zero variance or unexpressed) are excluded from the pool", {
  ref <- make_ref(n_genes = 30L, n_samp = 100L)
  ref[1, ] <- 5            # constant gene -> zero variance
  ref[2, ] <- 0            # unexpressed gene
  er <- build_empirical_reference(ref, p = 10, seed = 1)
  expect_equal(er$n_pool, 28L)
  expect_false("g1" %in% er$genes)
  expect_false("g2" %in% er$genes)
})

test_that("p beyond the usable pool samples with replacement and warns", {
  ref <- make_ref(n_genes = 8L, n_samp = 50L)
  expect_warning(
    er <- build_empirical_reference(ref, p = 20, seed = 1),
    "with replacement"
  )
  expect_equal(er$p, 20L)
  expect_true(anyDuplicated(er$genes) > 0)
})

test_that("invalid reference inputs are rejected", {
  expect_error(build_empirical_reference(1:10, p = 3), "numeric matrix")
  expect_error(build_empirical_reference(matrix(1, 5, 1), p = 2), "at least 2 samples")
})

test_that("low_var_features assigns the lowest-variance genes to the named features", {
  # gene g has sd ~ g, so variance increases with g
  n_genes <- 50L; n_samp <- 100L
  ref <- t(sapply(seq_len(n_genes), function(g) 10 + g * stats::qnorm(stats::ppoints(n_samp))))
  rownames(ref) <- paste0("g", seq_len(n_genes))

  lv_feats <- c(1L, 2L, 3L)
  er <- build_empirical_reference(ref, p = 20, seed = 1, low_var_features = lv_feats)

  v <- apply(er$quantiles, 2L, stats::var)
  expect_true(all(v[lv_feats] < min(v[-lv_feats])))   # named features hold the lowest-variance genes
})
