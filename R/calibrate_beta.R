#' Calibrate the global effect size to a target concordance
#'
#' Finds the `beta_scale` (the global relevance-strength knob of
#' [build_linear_predictor()]) at which the simulated data attain a target
#' Harrell concordance `target_C` on an uncensored pre-sample.
#'
#' @details
#' Under the proportional-hazards Bender-Weibull generator the uncensored
#' concordance depends only on the linear predictor: for a pair `(i, j)`,
#' `P(T_i < T_j | eta) = plogis(eta_i - eta_j)`, independent of the Weibull
#' shape (it cancels because the time map is monotone). Harrell's C is therefore
#' `E[plogis(|eta_i - eta_j|)]` over random pairs, which this function estimates
#' from `n_pairs` sampled pairs rather than all `O(n^2)` of them.
#'
#' Because `eta` scales linearly with `beta_scale`, the pre-sample and its
#' pairwise `|delta eta|` at unit scale are computed once, the bracket is seeded
#' with the heuristic `C ~ pnorm(sigma_eta / 2)` (`sigma_eta = sd(eta)`), and
#' `beta_scale` is then found by bisection. C is monotone increasing in
#' `beta_scale`, from `0.5` at `0` toward `1`.
#'
#' The pre-sample is drawn on the scenario's own marginal so the calibrated
#' concordance matches the data that will actually be simulated.
#'
#' @param Sigma,manifest The latent covariance and manifest from
#'   [make_sigma_blocks()].
#' @param target_C Target Harrell concordance, in `(0.5, 1)`. Default `0.70`.
#' @param marginals,negbin_mu,negbin_size,reference Passed to
#'   [sample_features()] to draw the pre-sample on the scenario marginal.
#' @param beta_int Relative interaction strength passed to
#'   [build_linear_predictor()]. Default `1`.
#' @param n_presample Pre-sample size. Default `4000`.
#' @param n_pairs Number of random pairs used to estimate C. Default `2e5`.
#' @param tol Tolerance on C for the bisection. Default `0.005`.
#' @param max_iter Maximum bisection iterations. Default `60`.
#' @param seed Integer or `NULL`. Seed for the reproducible pre-sample and pair
#'   draws.
#'
#' @return The calibrated `beta_scale` (a single numeric).
#'
#' @seealso [build_linear_predictor()], [sample_features()]
#' @export
calibrate_beta <- function(Sigma, manifest, target_C = 0.70,
                           marginals = c("gaussian", "negbin", "empirical"),
                           negbin_mu = 50, negbin_size = 2, reference = NULL,
                           beta_int = 1,
                           n_presample = 4000, n_pairs = 2e5,
                           tol = 0.005, max_iter = 60, seed = NULL) {
  marginals <- match.arg(marginals)
  stopifnot(
    "`target_C` must be a single number in (0.5, 1)" = length(target_C) == 1L && is.numeric(target_C) && target_C > 0.5 && target_C < 1,
    "`n_presample` must be at least 2"               = length(n_presample) == 1L && is.numeric(n_presample) && n_presample >= 2,
    "`n_pairs` must be at least 1"                    = length(n_pairs) == 1L && is.numeric(n_pairs) && n_pairs >= 1,
    "`tol` must be a single positive number"         = length(tol) == 1L && is.numeric(tol) && tol > 0
  )
  n_presample <- as.integer(n_presample)
  n_pairs <- as.integer(n_pairs)

  ## pre-sample once on the scenario marginal; collect pairwise |delta eta| at unit scale
  cal <- with_seed(seed, {
    X <- sample_features(n_presample, Sigma, marginals = marginals,
                         negbin_mu = negbin_mu, negbin_size = negbin_size,
                         reference = reference, seed = NULL, return_cor = FALSE)$X
    eta_unit <- build_linear_predictor(X, manifest, beta_scale = 1, beta_int = beta_int)
    i <- sample.int(n_presample, n_pairs, replace = TRUE)
    j <- sample.int(n_presample, n_pairs, replace = TRUE)
    same <- i == j
    j[same] <- (i[same] %% n_presample) + 1L
    list(d = abs(eta_unit[i] - eta_unit[j]), sd_eta = stats::sd(eta_unit))
  })
  if (!is.finite(cal$sd_eta) || cal$sd_eta == 0) {
    stop("eta has zero variance (no signal); cannot calibrate to C > 0.5.", call. = FALSE)
  }

  concordance_at <- function(s) mean(stats::plogis(s * cal$d))   # PH: P(T_i < T_j) = plogis(delta eta)
  bisect_scale(concordance_at, target_C, cal$sd_eta, tol = tol, max_iter = max_iter)
}

# internal: bisect the effect-scale `s` to hit `target_C`, shared by the PH and
# AFT calibrators. `concordance_at(s)` must increase from 0.5 (at s = 0) toward 1.
# The bracket starts from the PH heuristic C ~ pnorm(sigma_eta / 2) and doubles the
# upper bound until it exceeds the target, so it also brackets the AFT arm.
bisect_scale <- function(concordance_at, target_C, sd_eta, tol = 0.005, max_iter = 60) {
  s_hi <- 2 * stats::qnorm(target_C) / sd_eta
  if (!is.finite(s_hi) || s_hi <= 0) s_hi <- 1
  expand <- 0L
  while (concordance_at(s_hi) < target_C && expand < 60L) {
    s_hi <- s_hi * 2
    expand <- expand + 1L
  }
  lo <- 0
  hi <- s_hi
  mid <- hi
  for (it in seq_len(max_iter)) {
    mid <- (lo + hi) / 2
    cm <- concordance_at(mid)
    if (abs(cm - target_C) < tol) break
    if (cm < target_C) lo <- mid else hi <- mid
  }
  mid
}

# internal: empirical effect-size calibration for the log-logistic AFT arm. The
# PH plogis formula in calibrate_beta() is exact only under proportional hazards,
# so the AFT arm would otherwise overshoot target_C. The AFT log-time is
# logT = logT0 - eta with eta = s * eta_unit and logT0 = qlogis(U) / shape; fixing
# the baseline noise logT0 once makes the pairwise concordance a deterministic,
# (effectively monotone) function of s, calibrated by the shared bisection. The
# pre-sample is drawn on the scenario marginal, like calibrate_beta().
calibrate_beta_aft <- function(Sigma, manifest, target_C = 0.70,
                               marginals = "gaussian", negbin_mu = 50, negbin_size = 2,
                               reference = NULL, beta_int = 1, aft_shape = 3,
                               n_presample = 4000, n_pairs = 2e5,
                               tol = 0.005, max_iter = 60, seed = NULL) {
  n_presample <- as.integer(n_presample)
  n_pairs <- as.integer(n_pairs)

  cal <- with_seed(seed, {
    X <- sample_features(n_presample, Sigma, marginals = marginals,
                         negbin_mu = negbin_mu, negbin_size = negbin_size,
                         reference = reference, seed = NULL, return_cor = FALSE)$X
    eta_unit <- build_linear_predictor(X, manifest, beta_scale = 1, beta_int = beta_int)
    logT0 <- stats::qlogis(stats::runif(n_presample)) / aft_shape   # fixed baseline AFT log-times
    i <- sample.int(n_presample, n_pairs, replace = TRUE)
    j <- sample.int(n_presample, n_pairs, replace = TRUE)
    same <- i == j
    j[same] <- (i[same] %% n_presample) + 1L
    list(d_eta = eta_unit[i] - eta_unit[j], d_logT0 = logT0[i] - logT0[j],
         sd_eta = stats::sd(eta_unit))
  })
  if (!is.finite(cal$sd_eta) || cal$sd_eta == 0) {
    stop("eta has zero variance (no signal); cannot calibrate to C > 0.5.", call. = FALSE)
  }

  ## pairwise AFT concordance at scale s (baseline noise fixed):
  ## eta_i - eta_j = s * d_eta ; logT_i - logT_j = d_logT0 - s * d_eta (time is monotone in logT)
  concordance_at <- function(s) {
    de <- s * cal$d_eta
    dt <- cal$d_logT0 - s * cal$d_eta
    comparable <- de != 0 & dt != 0
    if (!any(comparable)) return(0.5)
    concordant <- (de > 0 & dt < 0) | (de < 0 & dt > 0)   # higher eta, shorter time
    mean(concordant[comparable])
  }
  bisect_scale(concordance_at, target_C, cal$sd_eta, tol = tol, max_iter = max_iter)
}
