#' Default simulation parameters
#'
#' Returns the default simulation scenario as a single config list, with any
#' overrides applied. Pass this (optionally overridden) to [simulate_dataset()].
#'
#' @details
#' Fields: structure (`p`, `k_drivers`, `cluster_size`, `rho`,
#' `n_distractor_subblocks`, `distractor_size`, `distractor_rho`,
#' `interaction`); sampling (`n`, `marginals`, `negbin_mu`, `negbin_size`,
#' `reference`); effect (`beta_int`, `target_C`, `driver_var_scale`); event times (`time_model` one
#' of `"weibull_ph"`/`"loglogistic_aft"`, `weibull_shape`, `weibull_lambda`,
#' `aft_shape`); censoring (`censor_fraction`, `censor_mechanism`,
#' `fraction_model`, `gamma`, `disjoint_cols`); reproducibility (`seed`); and
#' calibration controls (`calib_n_presample`, `calib_n_pairs`).
#'
#' @param ... Named overrides, e.g. `default_params(n = 400, censor_fraction = 0.7)`.
#'
#' @return A named list of simulation parameters.
#' @seealso [simulate_dataset()]
#'
#' @examples
#' # override only what differs from the defaults
#' params <- default_params(p = 60, n = 200, censor_fraction = 0.5, seed = 1)
#' params[c("p", "n", "censor_fraction")]
#' @export
default_params <- function(...) {
  p <- list(
    p = 5000L, k_drivers = 8L, cluster_size = 5L, rho = 0.8,
    n_distractor_subblocks = 20L, distractor_size = 10L, distractor_rho = 0.7,
    interaction = TRUE,
    n = 200L,
    marginals = "gaussian", negbin_mu = 50, negbin_size = 2, reference = NULL,
    beta_int = 1, target_C = 0.70, driver_var_scale = 1,
    time_model = "weibull_ph", weibull_shape = 1.5, weibull_lambda = log(2), aft_shape = 3,
    censor_fraction = 0.5, censor_mechanism = "independent", fraction_model = "exponential",
    gamma = 0.7, disjoint_cols = NULL,
    seed = NULL, calib_n_presample = 4000L, calib_n_pairs = 2e5
  )
  ov <- list(...)
  if (length(ov)) p[names(ov)] <- ov
  p
}

# internal: switch to an L'Ecuyer-CMRG stream seeded by `seed`, returning a thunk
# that restores the previous RNG kind and state. With seed = NULL nothing is
# changed (the current, possibly sweep-managed, stream is used).
setup_lecuyer <- function(seed) {
  if (is.null(seed)) return(function() invisible(NULL))
  old_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  function() {
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
    if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }
}

# internal: exact uncensored Harrell concordance of the event times (arm-agnostic:
# correct for both the PH and the AFT arms, unlike the PH-only plogis formula).
empirical_concordance <- function(eta, time) {
  de <- outer(eta, eta, "-")
  dt <- outer(time, time, "-")
  up <- upper.tri(de)
  de <- de[up]
  dt <- dt[up]
  comparable <- de != 0 & dt != 0
  concordant <- (de > 0 & dt < 0) | (de < 0 & dt > 0)   # higher eta, shorter time
  mean(concordant[comparable])
}

# internal: mean realised intra-cluster correlation across the driver clusters
realized_intracluster_cor <- function(X, manifest) {
  ids <- unique(manifest$cluster_id[manifest$role == "driver"])
  vals <- vapply(ids, function(cid) {
    idx <- which(manifest$cluster_id == cid & manifest$role %in% c("driver", "redundant"))
    if (length(idx) < 2L) return(NA_real_)
    cc <- stats::cor(X[, idx])
    mean(cc[upper.tri(cc)])
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Simulate one right-censored survival replicate with known ground truth
#'
#' Composes the simulator's building blocks into a single replicate:
#' block latent covariance, copula feature sampling, the hazard linear
#' predictor, effect-size calibration, event times, and censoring. Returns the
#' observed data, the truth manifest, and realised diagnostics.
#'
#' @details
#' Random draws are seeded once with an L'Ecuyer-CMRG stream from `params$seed`
#' (restored on exit); set `params$seed = NULL` to draw from the current stream,
#' as the sweep layer does when it threads per-replicate streams. Effect size is
#' calibrated to `params$target_C` on the scenario marginal (the PH arm
#' analytically, the AFT arm by empirical bisection on its own concordance, since
#' the analytic formula holds only under proportional hazards), and the censoring
#' parameter to `params$censor_fraction`. When `params$driver_var_scale < 1` the
#' driver-cluster columns of the returned `X` are shrunk to a low marginal
#' variance (rank-preserving), the low-variance-signal scenario in which a
#' variance filter fails but FCBF does not.
#'
#' @param params A config list, normally from [default_params()].
#'
#' @return A list with `X` (`n` by `p` features), `time`, `status`, `manifest`,
#'   and `realized` (a list with `beta_scale`, `censoring_fraction`, `n_events`,
#'   the uncensored `concordance`, and the realised `intracluster_cor`).
#'
#' @seealso [default_params()]
#'
#' @examples
#' sim <- simulate_dataset(default_params(
#'   p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'   n_distractor_subblocks = 2L, distractor_size = 4L,
#'   calib_n_presample = 500, seed = 1))
#' dim(sim$X)
#' table(sim$manifest$role)   # the known truth, one row per feature
#' mean(sim$status)           # realised event fraction
#' @export
simulate_dataset <- function(params = default_params()) {
  stopifnot("`params` must be a list" = is.list(params))
  pp <- params
  restore <- setup_lecuyer(pp$seed)
  on.exit(restore(), add = TRUE)

  blocks <- make_sigma_blocks(pp$p, pp$k_drivers, pp$cluster_size, pp$rho,
                              pp$n_distractor_subblocks, pp$distractor_size, pp$distractor_rho,
                              interaction = pp$interaction)

  ## Effect-size calibration is arm-specific: the PH arm uses the exact analytic
  ## plogis concordance, the AFT arm an empirical bisection on its own generator
  ## (the plogis formula only holds under proportional hazards).
  beta_scale <- if (identical(pp$time_model, "loglogistic_aft")) {
    calibrate_beta_aft(blocks$Sigma, blocks$manifest, target_C = pp$target_C,
                       marginals = pp$marginals, negbin_mu = pp$negbin_mu,
                       negbin_size = pp$negbin_size, reference = pp$reference,
                       beta_int = pp$beta_int, aft_shape = pp$aft_shape,
                       n_presample = pp$calib_n_presample, n_pairs = pp$calib_n_pairs, seed = NULL)
  } else {
    calibrate_beta(blocks$Sigma, blocks$manifest, target_C = pp$target_C,
                   marginals = pp$marginals, negbin_mu = pp$negbin_mu,
                   negbin_size = pp$negbin_size, reference = pp$reference,
                   beta_int = pp$beta_int, n_presample = pp$calib_n_presample,
                   n_pairs = pp$calib_n_pairs, seed = NULL)
  }

  X <- sample_features(pp$n, blocks$Sigma, marginals = pp$marginals,
                       negbin_mu = pp$negbin_mu, negbin_size = pp$negbin_size,
                       reference = pp$reference, seed = NULL, return_cor = FALSE)$X

  eta <- build_linear_predictor(X, blocks$manifest, beta_scale = beta_scale, beta_int = pp$beta_int)

  Tt <- switch(pp$time_model,
    weibull_ph = bender_weibull_time(eta, shape_a = pp$weibull_shape,
                                     scale_lambda = pp$weibull_lambda, seed = NULL),
    loglogistic_aft = aft_loglogistic_time(eta, shape = pp$aft_shape, seed = NULL),
    stop("unknown `time_model`: ", pp$time_model, call. = FALSE)
  )

  cens_param <- calibrate_censoring(eta, X, Tt, target_fraction = pp$censor_fraction,
                                    mechanism = pp$censor_mechanism, fraction_model = pp$fraction_model,
                                    manifest = blocks$manifest, gamma = pp$gamma,
                                    disjoint_cols = pp$disjoint_cols, seed = NULL)
  C <- make_censoring(eta, X, cens_param, mechanism = pp$censor_mechanism,
                      fraction_model = pp$fraction_model, manifest = blocks$manifest,
                      gamma = pp$gamma, disjoint_cols = pp$disjoint_cols, seed = NULL)

  time <- pmin(Tt, C)
  status <- as.integer(Tt <= C)

  ## low-variance-signal scenario: shrink the driver-cluster marginal variance.
  ## Scaling is rank-preserving, so the FCBF relevance/redundancy ordering is
  ## unchanged, but a variance filter is misled into ranking the signal low.
  if (pp$driver_var_scale != 1) {
    sig_cols <- which(blocks$manifest$role %in% c("driver", "redundant"))
    X[, sig_cols] <- X[, sig_cols] * pp$driver_var_scale
  }

  realized <- list(
    beta_scale = beta_scale,
    censoring_fraction = mean(status == 0L),
    n_events = sum(status),
    concordance = empirical_concordance(eta, Tt),
    intracluster_cor = realized_intracluster_cor(X, blocks$manifest)
  )

  list(X = X, time = time, status = status, manifest = blocks$manifest, realized = realized)
}
