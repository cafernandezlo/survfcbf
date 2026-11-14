#' Leakage-free FCBF feature selection for right-censored survival data
#'
#' The survival adaptation of FCBF: on a single (training)
#' fold it transforms the censored outcome to a pseudo-target, discretizes the
#' target and every feature, sets the relevance threshold, and runs the FCBF
#' predominance selection. Everything learned from data (the null Cox model, the
#' per-feature cutpoints, `delta`) comes only from the supplied fold, and the
#' learned cutpoints are returned so a held-out fold can be discretized
#' identically, with no leakage.
#'
#' @details
#' Pipeline: [survival_target()] (martingale residuals of a null Cox model) ->
#' [discretize_quantile()] of the target and each feature into `q` equal-frequency
#' bins -> `delta` from [delta_permutation()] (unless a fixed `delta` is supplied,
#' the inner-CV-grid alternative) -> [fcbf_select()]. Only the permutation null
#' uses random numbers, so `seed` makes the whole call reproducible; with a fixed
#' `delta` it is deterministic.
#'
#' Supplying `top` switches to the explicit-cardinality mode: the relevance/redundancy
#' machinery is unchanged (SU relevance, SU Markov-blanket pruning), but the set size
#' is capped at `top` predominant features instead of being decided by `delta`, and the
#' permutation null is skipped (so the call is deterministic without a `seed`). This is
#' the mode for the low-event / high-dimensional genomic regime, where the
#' permutation-null `delta` inflates above the signal and self-selection collapses
#' (`q = 2` is the validated choice there). It keeps the FCBF identity; only the
#' cardinality control changes.
#'
#' Coarse discretisation makes many features share the same relevance (at `q = 2`
#' SU takes only a few dozen distinct values over thousands of genes). By default
#' those ties are broken by the absolute Spearman correlation of each continuous
#' feature with the pseudo-target, computed on the same training data, so the
#' selection does not depend on the column order of `X`. `tiebreak = "index"`
#' restores the classic column-order rule of the original implementation.
#'
#' @section Censoring in the top target bin:
#' An event has a martingale residual `1 - H0(t)` and a censored patient `-H0(t)`, so
#' the censored patients with the shortest follow-up rank just below the events. When
#' the events are fewer than the patients in the top bin of the discretized target, as
#' with `q = 2` and more than half of the patients censored, those short-follow-up
#' censored patients fill the rest of the high-risk class. The relevance of a gene can
#' then partly reflect follow-up length rather than risk.
#'
#' When more than half of the patients in the top bin are censored,
#' `survfcbf_select()` raises a warning of class `survfcbf_followup_warning` that
#' reports the censored share. It is computed on the data passed to the function (the
#' training fold), in top-K and threshold mode alike and for any `q`. The warning only
#' informs and the selection is the same with or without it. It is worth checking the
#' result, for example by rerunning with a larger `q`. Because the warning has its own
#' class it can be silenced alone with [withCallingHandlers()] and
#' `invokeRestart("muffleWarning")`, as in the examples.
#'
#' The check covers the martingale target with SU relevance, the only setting where
#' censored patients enter the discretized target that defines relevance. Under
#' `target_method = "ipcw_binned"` censored patients carry zero weight, and under
#' `relevance = "cox_score"` relevance is a concordance index that handles censoring
#' itself, so neither raises the warning.
#'
#' @param time,status Right-censored outcome for the training fold.
#' @param X Numeric `n` by `p` feature matrix for the training fold.
#' @param q Number of discretization bins. `NULL` (default) uses Bommert's rule.
#' @param B Permutations for the `delta` null. Default `200`.
#' @param target_method Target transform passed to [survival_target()]:
#'   `"martingale"` (default) or `"ipcw_binned"` (the 5.5 ablation, which uses the
#'   observed time as the target with IPCW-weighted SU throughout).
#' @param relevance Relevance measure: `"su"` (symmetrical uncertainty, default)
#'   or `"cox_score"` (rescaled univariate Cox concordance `2 |C - 0.5|`, the 5.5
#'   ablation; redundancy still uses SU, and an explicit `delta` is required).
#' @param delta Fixed relevance threshold; `NULL` (default) estimates it by the
#'   permutation null. Under `top` it acts only as an optional relevance floor.
#' @param seed Integer or `NULL`. Seed for the permutation null.
#' @param top Optional positive integer. When supplied, selects the `top` most
#'   relevant predominant features (explicit-cardinality mode) instead of letting
#'   `delta` determine the set size, and skips the permutation null. `NULL`
#'   (default) is the classic delta-driven selection.
#' @param tiebreak How features with equal relevance are ordered: `"spearman"`
#'   (default), by the absolute Spearman correlation of the continuous feature with
#'   the pseudo-target, or `"index"`, by column position. Only `"spearman"` makes
#'   the selection invariant to the column order of `X`.
#'
#' @return A list with `selected` (predominant feature column indices),
#'   `cutpoints` (length-`p` list of per-feature break vectors for re-discretizing
#'   held-out data), `delta` (the threshold used, or `NULL` under `top` when no
#'   floor was supplied), `target` (the training pseudo-target, `NULL` under
#'   `cox_score`), `relevance` (per-feature relevance: `SU(Fi, y)`, or the rescaled
#'   Cox score under `relevance = "cox_score"`) and `trace` (the FCBF
#'   Markov-blanket trace; its `relevance_su` slot holds whichever relevance was
#'   used, so it is the Cox score rather than SU under that ablation).
#'
#' @seealso [survival_target()], [discretize_quantile()], [delta_permutation()],
#'   [fcbf_select()]
#'
#' @examples
#' sim <- simulate_dataset(default_params(
#'   p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'   n_distractor_subblocks = 2L, distractor_size = 4L,
#'   calib_n_presample = 500, seed = 1))
#'
#' # top-K mode, the recommended one: deterministic, no permutation null
#' sel <- survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 5)
#' sel$selected
#'
#' # the same genes whatever the column order
#' perm <- sample(ncol(sim$X))
#' sel2 <- survfcbf_select(sim$time, sim$status, sim$X[, perm], q = 2, top = 5)
#' setequal(perm[sel2$selected], sel$selected)
#'
#' # recovery of the planted driver clusters
#' score_selection_clusterwise(sel$selected, sim$manifest)$sensitivity
#'
#' # silence only the follow-up warning, any other warning still shows
#' quiet <- withCallingHandlers(
#'   survfcbf_select(sim$time, sim$status, sim$X, q = 2, top = 5),
#'   survfcbf_followup_warning = function(w) invokeRestart("muffleWarning"))
#' identical(quiet$selected, sel$selected)
#' @export
survfcbf_select <- function(time, status, X, q = NULL, B = 200L,
                            target_method = c("martingale", "ipcw_binned"),
                            relevance = c("su", "cox_score"),
                            delta = NULL, seed = NULL, top = NULL,
                            tiebreak = c("spearman", "index")) {
  relevance <- match.arg(relevance)
  target_method <- match.arg(target_method)
  tiebreak <- match.arg(tiebreak)
  stopifnot(
    "`X` must be a numeric matrix"                 = is.matrix(X) && is.numeric(X),
    "`X`, `time`, `status` must share nrow/length" = nrow(X) == length(time) && length(time) == length(status),
    "`top` must be NULL or a positive integer"     = is.null(top) || (length(top) == 1L && is.numeric(top) && !is.na(top) && top >= 1)
  )
  if (!is.null(top)) top <- as.integer(top)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(q)) {
    q <- max(min(floor(n / 3), 10L), 2L)
  }

  ## discretise the features (cut-points learned on train), used by both relevance legs
  feat <- lapply(seq_len(p), function(j) discretize_quantile(X[, j], q = q))
  X_disc <- vapply(feat, function(f) f$binned, integer(n))
  cutpoints <- lapply(feat, function(f) f$breaks)
  followup <- NULL   # censored share of the top target bin, checked for the martingale target only

  if (relevance == "su") {
    ## pseudo-target -> discretise -> delta (permutation null unless fixed) -> FCBF.
    ## The ipcw_binned target uses the observed time with IPCW-weighted SU throughout.
    ym <- survival_target(time, status, method = target_method)
    w <- if (target_method == "ipcw_binned") ipcw_weights(time, status) else NULL
    y_disc <- discretize_quantile(ym, q = q)$binned
    if (target_method == "martingale") followup <- top_bin_censoring(y_disc, status)
    ## the permutation-null delta is only the cardinality control for delta mode;
    ## under `top` the explicit cap selects the cardinality and the null is skipped.
    if (is.null(top) && is.null(delta)) {
      delta <- delta_permutation(X_disc, y_disc, B = B, weights = w, seed = seed)
    }
    tb <- if (tiebreak == "spearman") spearman_tiebreak(X, ym) else NULL
    fc <- fcbf_select(X_disc, y_disc, delta = if (is.null(delta)) 0 else delta, weights = w, top = top,
                      tiebreak = tb)
    rel <- fc$trace$relevance_su
    target <- ym
  } else {
    ## Cox-score ablation: rescaled univariate concordance as the relevance leg,
    ## SU still for redundancy. Its permutation null is impractical at scale, so the
    ## cardinality must come from an explicit `delta` (fixed-grid) or a `top` cap.
    if (is.null(top) && is.null(delta)) {
      stop("relevance = 'cox_score' needs an explicit `delta` or a `top` cardinality (its permutation null is impractical at scale; use the fixed-grid alternative or top-K).",
           call. = FALSE)
    }
    rel <- cox_score_relevance(time, status, X)
    tb <- if (tiebreak == "spearman") {
      spearman_tiebreak(X, survival_target(time, status, method = "martingale"))
    } else NULL
    fc <- fcbf_select(X_disc, delta = if (is.null(delta)) 0 else delta, relevance = rel, top = top,
                      tiebreak = tb)
    target <- NULL
  }

  ## raised after the selection is complete and from data already computed, so it
  ## cannot change the result
  if (!is.null(followup) && followup$prop > 0.5) followup_warning(followup, q)

  list(
    selected = fc$selected,
    cutpoints = cutpoints,
    delta = delta,
    target = target,
    relevance = rel,
    trace = fc$trace
  )
}

# internal: censored share of the top bin of the discretized target. The top bin is
# the highest realised code, so it stays right when duplicate cutpoints collapse bins.
top_bin_censoring <- function(y_disc, status) {
  top <- y_disc == max(y_disc)
  n_top <- sum(top)
  n_censored <- sum(top & status == 0)
  list(n_top = n_top, n_censored = n_censored, prop = n_censored / n_top)
}

# internal: the follow-up warning, with its own class so it can be muffled alone.
# The counts travel in the condition object as well as in the message.
followup_warning <- function(followup, q) {
  msg <- sprintf(paste(
    "%.1f%% of the patients in the top bin of the discretized survival target are censored (%d of %d, q = %d).",
    "Feature relevance may partly reflect follow-up length rather than risk.",
    "Check the selection, for example with a larger q."),
    100 * followup$prop, as.integer(followup$n_censored), as.integer(followup$n_top), as.integer(q))
  warning(warningCondition(msg, prop = followup$prop, n_top = followup$n_top,
                           n_censored = followup$n_censored,
                           class = "survfcbf_followup_warning", call = NULL))
}

# internal: secondary ordering key for features tied on relevance, the absolute
# Spearman correlation of each continuous feature with the pseudo-target. A constant
# feature has no correlation and scores 0.
spearman_tiebreak <- function(X, target) {
  r <- suppressWarnings(stats::cor(X, target, method = "spearman"))
  r <- abs(as.numeric(r))
  r[!is.finite(r)] <- 0
  r
}

# internal: IPCW weights for the ipcw_binned target. Events are weighted by
# 1 / G_hat(time) with G_hat the Kaplan-Meier of the censoring distribution;
# censored subjects get weight 0 (the standard complete-case IPCW scheme).
ipcw_weights <- function(time, status) {
  if (sum(status == 1) == 0L) {
    stop("ipcw_binned needs at least one event (all observations are censored).", call. = FALSE)
  }
  cfit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)   # KM of the censoring distribution
  # Left-continuous censoring survival G(t-): its value just before t, i.e. at the
  # last censoring knot strictly below t (1 before the first). Using G(t-) instead
  # of G(t) avoids over-weighting events that tie a censoring time.
  j <- findInterval(time, cfit$time, left.open = TRUE)             # number of knots strictly < t
  G <- ifelse(j == 0L, 1, cfit$surv[pmax(j, 1L)])
  ifelse(status == 1, 1 / pmax(G, 1e-6), 0)
}

# internal: univariate Cox-score relevance, rescaled to [0, 1] (SU-comparable).
# 2 * |C - 0.5| of coxph(Surv ~ feature) per feature; a constant feature scores 0.
cox_score_relevance <- function(time, status, X) {
  surv <- survival::Surv(time, status)
  vapply(seq_len(ncol(X)), function(j) {
    xj <- X[, j]
    if (stats::sd(xj) == 0) return(0)
    cc <- survival::concordance(surv ~ xj)$concordance
    2 * abs(cc - 0.5)
  }, numeric(1))
}
