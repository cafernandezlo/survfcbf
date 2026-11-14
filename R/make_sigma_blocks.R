#' Build a block-structured latent covariance matrix with known ground truth
#'
#' Constructs the latent Gaussian covariance `Sigma` (`p` by `p`) used by the
#' survival simulator together with a `manifest` describing
#' the role, cluster membership and hazard coefficient of every column. The
#' structure embeds *designed* redundancy so the approximate Markov blanket /
#' predominance logic of FCBF has something to act on, and so that recovery can
#' be scored at cluster granularity rather than by exact feature index.
#'
#' @details
#' `Sigma` is block diagonal across five kinds of column, all with unit latent
#' variance (so `Sigma` is a correlation matrix and relevance is not confounded
#' with marginal variance):
#'
#' 1. **Drivers.** `k_drivers` features, each carrying a non-zero `beta`. They
#'    are mutually near-independent because each driver is the first member of
#'    its own redundant cluster block and the blocks are independent.
#' 2. **Redundant clusters.** One compound-symmetry block per driver,
#'    `Sigma_block = (1 - rho) I + rho 11'`, spanning the driver *and* its
#'    `cluster_size` redundants. A cluster block therefore has
#'    `cluster_size + 1` columns. Only the driver carries the effect; the
#'    redundants have conditional `beta = 0`, so FCBF should keep one cluster
#'    member and drop the rest.
#' 3. **Distractors.** `n_distractor_subblocks` equicorrelated subblocks of
#'    `distractor_size` columns at correlation `distractor_rho`, all with
#'    `beta = 0`. Correlated but irrelevant: they measure false-discovery
#'    control under correlation.
#' 4. **Noise.** Pure i.i.d. columns filling the remainder up to `p`, forcing
#'    `p >> n`.
#' 5. **Interaction.** If `interaction`, two independent columns with marginal
#'    `beta = 0` whose product drives the hazard (the product term is set later,
#'    in `build_linear_predictor`). FCBF, being univariate in its relevance leg,
#'    is expected to miss these; they are included to exhibit the limitation
#'    honestly.
#'
#' Columns are laid out in the order clusters, distractors, interaction, noise.
#' Row `i` of `manifest` corresponds to column `i` of `Sigma`.
#'
#' Driver coefficients carry mixed signs (balanced alternating unit magnitudes
#' `+1, -1, +1, ...`); the global effect size is applied later through
#' `beta_scale` in `build_linear_predictor`. Every other role has `beta = 0`.
#'
#' Each compound-symmetry block of size `m` and correlation `r` has eigenvalues
#' `1 + (m - 1) r` (once) and `1 - r` (`m - 1` times), so it is positive definite
#' iff `1 - r > 0` and `1 + (m - 1) r > 0`. The function stops with an error if
#' `rho` or `distractor_rho` violate this, guaranteeing a positive-definite
#' `Sigma`.
#'
#' This function is deterministic given its arguments: no random numbers are
#' drawn here (sampling happens in `sample_features`).
#'
#' @param p Integer. Total number of features (columns of `Sigma`). Must be at
#'   least the number of structured columns,
#'   `k_drivers * (cluster_size + 1) + n_distractor_subblocks * distractor_size`
#'   plus `2` if `interaction` is `TRUE`; the remainder becomes noise.
#' @param k_drivers Integer. Number of driver features, one per redundant
#'   cluster.
#' @param cluster_size Integer. Number of redundants per driver. The
#'   corresponding cluster block has `cluster_size + 1` columns (driver
#'   included).
#' @param rho Numeric. Intra-cluster latent correlation for the driver and
#'   redundant compound-symmetry blocks.
#' @param n_distractor_subblocks Integer. Number of equicorrelated distractor
#'   subblocks.
#' @param distractor_size Integer. Number of columns per distractor subblock.
#' @param distractor_rho Numeric. Equicorrelation within each distractor
#'   subblock.
#' @param interaction Logical. If `TRUE` (default) add an independent pair of
#'   interaction columns with marginal `beta = 0`.
#'
#' @return A named list with two elements:
#'
#'   * `Sigma`: the `p` by `p` symmetric positive-definite latent correlation
#'     matrix (unit diagonal), block diagonal across all groups.
#'   * `manifest`: a `data.frame` with one row per column of `Sigma`, in column
#'     order, and three columns: `role` (one of `"driver"`, `"redundant"`,
#'     `"distractor"`, `"noise"`, `"interaction"`), `cluster_id` (integer group
#'     id shared within a driver cluster and within a distractor subblock, `NA`
#'     for noise), and `beta` (the latent hazard coefficient, non-zero only for
#'     drivers).
#'
#' @examples
#' blocks <- make_sigma_blocks(
#'   p = 60, k_drivers = 3, cluster_size = 4, rho = 0.8,
#'   n_distractor_subblocks = 2, distractor_size = 5, distractor_rho = 0.7
#' )
#' dim(blocks$Sigma)
#' table(blocks$manifest$role)
#'
#' @export
make_sigma_blocks <- function(p,
                              k_drivers,
                              cluster_size,
                              rho,
                              n_distractor_subblocks,
                              distractor_size,
                              distractor_rho,
                              interaction = TRUE) {
  ## ---- argument checks ------------------------------------------------------
  is_count <- function(x) {
    length(x) == 1L && is.numeric(x) && is.finite(x) && x == as.integer(x)
  }
  is_num1 <- function(x) length(x) == 1L && is.numeric(x) && is.finite(x)
  stopifnot(
    "`p` must be a single positive integer"                   = is_count(p) && p >= 1,
    "`k_drivers` must be a single positive integer"           = is_count(k_drivers) && k_drivers >= 1,
    "`cluster_size` must be a single non-negative integer"    = is_count(cluster_size) && cluster_size >= 0,
    "`n_distractor_subblocks` must be a non-negative integer" = is_count(n_distractor_subblocks) && n_distractor_subblocks >= 0,
    "`distractor_size` must be a single positive integer"     = is_count(distractor_size) && distractor_size >= 1,
    "`rho` must be a single finite number"                    = is_num1(rho),
    "`distractor_rho` must be a single finite number"         = is_num1(distractor_rho),
    "`interaction` must be a single non-NA logical"           = length(interaction) == 1L && is.logical(interaction) && !is.na(interaction)
  )

  p                      <- as.integer(p)
  k_drivers              <- as.integer(k_drivers)
  cluster_size           <- as.integer(cluster_size)
  n_distractor_subblocks <- as.integer(n_distractor_subblocks)
  distractor_size        <- as.integer(distractor_size)

  ## ---- positive-definiteness guard for compound-symmetry blocks -------------
  ## A size-m equicorrelation block (1 - r) I + r 11' is PD iff 1 - r > 0 and
  ## 1 + (m - 1) r > 0. A 1x1 block is just [1] and always PD.
  cs_is_pd <- function(m, r) {
    if (m <= 1L) return(TRUE)
    (1 - r) > 0 && (1 + (m - 1L) * r) > 0
  }
  cluster_block_size <- cluster_size + 1L            # driver + its redundants
  if (!cs_is_pd(cluster_block_size, rho)) {
    stop(sprintf(
      "`rho` = %g does not keep a %d-column cluster block positive definite; require -1/%d < rho < 1.",
      rho, cluster_block_size, cluster_block_size - 1L
    ), call. = FALSE)
  }
  if (!cs_is_pd(distractor_size, distractor_rho)) {
    stop(sprintf(
      "`distractor_rho` = %g does not keep a %d-column distractor block positive definite; require -1/%d < distractor_rho < 1.",
      distractor_rho, distractor_size, max(distractor_size - 1L, 1L)
    ), call. = FALSE)
  }

  ## ---- block sizes and layout -----------------------------------------------
  n_cluster_cols     <- k_drivers * cluster_block_size
  n_distractor_cols  <- n_distractor_subblocks * distractor_size
  n_interaction_cols <- if (interaction) 2L else 0L
  n_structured       <- n_cluster_cols + n_distractor_cols + n_interaction_cols
  n_noise            <- p - n_structured
  if (n_noise < 0L) {
    stop(sprintf(
      "`p` = %d is too small: structured columns already require %d (clusters %d + distractors %d + interaction %d).",
      p, n_structured, n_cluster_cols, n_distractor_cols, n_interaction_cols
    ), call. = FALSE)
  }

  ## ---- compound-symmetry block builder --------------------------------------
  cs_block <- function(m, r) {
    B <- matrix(r, m, m)
    diag(B) <- 1
    B
  }

  ## ---- assemble Sigma and the manifest --------------------------------------
  Sigma <- diag(nrow = p)                 # unit variances; noise + interaction stay i.i.d.

  role       <- character(p)
  cluster_id <- rep(NA_integer_, p)
  beta       <- numeric(p)                # 0 everywhere unless set below

  ## driver coefficients: balanced alternating unit magnitudes (mixed signs)
  driver_beta <- rep(c(1, -1), length.out = k_drivers)

  pos <- 0L                               # running column offset

  ## (1) + (2) driver/redundant clusters
  for (cl in seq_len(k_drivers)) {
    idx <- pos + seq_len(cluster_block_size)
    Sigma[idx, idx] <- cs_block(cluster_block_size, rho)
    role[idx]       <- c("driver", rep("redundant", cluster_size))
    cluster_id[idx] <- cl
    beta[idx[1L]]   <- driver_beta[cl]    # only the driver carries the effect
    pos <- pos + cluster_block_size
  }

  ## (3) distractor subblocks
  for (s in seq_len(n_distractor_subblocks)) {
    idx <- pos + seq_len(distractor_size)
    Sigma[idx, idx] <- cs_block(distractor_size, distractor_rho)
    role[idx]       <- "distractor"
    cluster_id[idx] <- k_drivers + s
    pos <- pos + distractor_size
  }

  ## (5) interaction pair (independent; Sigma block already identity)
  if (interaction) {
    idx <- pos + seq_len(2L)
    role[idx]       <- "interaction"
    cluster_id[idx] <- k_drivers + n_distractor_subblocks + 1L
    pos <- pos + 2L
  }

  ## (4) noise (Sigma block already identity; cluster_id stays NA)
  if (n_noise > 0L) {
    idx <- pos + seq_len(n_noise)
    role[idx] <- "noise"
    pos <- pos + n_noise
  }

  manifest <- data.frame(
    role       = role,
    cluster_id = cluster_id,
    beta       = beta,
    stringsAsFactors = FALSE
  )

  list(Sigma = Sigma, manifest = manifest)
}
