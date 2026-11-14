#' Score a selected feature set against the cluster-level ground truth
#'
#' Evaluates a feature selection at *cluster* granularity,
#' the truth scale that matters when redundant features are interchangeable: a
#' relevant cluster counts as recovered as soon as any one of its members is
#' selected, distractor and noise selections are false discoveries, and
#' selecting more than one member of a relevant cluster is charged to parsimony,
#' not sensitivity.
#'
#' @details
#' Relevant clusters are the driver clusters (one per driver, the driver plus
#' its redundants sharing a `cluster_id`). Interaction features are a separate
#' category: FCBF's univariate relevance leg is expected to miss them, so they
#' count neither as true positives nor as false discoveries (they only enter the
#' false-discovery denominator as ordinary selections).
#'
#' Metrics:
#' * `sensitivity`: recovered relevant clusters divided by the number of
#'   relevant clusters.
#' * `fdr`: distractor-or-noise selections divided by the total number of
#'   selections (`0` when nothing is selected).
#' * `f1`: cluster-collapsed F1 with `precision = TP / (TP + FP)` and
#'   `recall = sensitivity`, where `TP` is recovered relevant clusters and `FP`
#'   is distractor-or-noise selections.
#' * `exact_recovery`: `TRUE` when every relevant cluster is recovered and no
#'   distractor or noise feature is selected (parsimony aside).
#' * `parsimony_penalty`: summed over relevant clusters, the redundant
#'   within-cluster selections beyond the first.
#' * `representative_per_cluster`: for each recovered cluster, the selected
#'   member indices (the driver first), for representative-stability tracking.
#'
#' @param selected Selected features, either integer column indices into
#'   `manifest` or a logical vector of length `nrow(manifest)`.
#' @param manifest The [make_sigma_blocks()] manifest, with `role` and
#'   `cluster_id` columns.
#'
#' @return A named list `sensitivity`, `fdr`, `f1`, `exact_recovery`,
#'   `parsimony_penalty`, `representative_per_cluster`.
#'
#' @seealso [make_sigma_blocks()]
#'
#' @examples
#' sim <- simulate_dataset(default_params(
#'   p = 60, n = 200, k_drivers = 3L, cluster_size = 3L,
#'   n_distractor_subblocks = 2L, distractor_size = 4L,
#'   calib_n_presample = 500, seed = 1))
#' drivers <- which(sim$manifest$role == "driver")
#'
#' # one member per driver cluster and nothing else: exact recovery
#' score_selection_clusterwise(drivers, sim$manifest)$exact_recovery
#'
#' # any member of a redundant cluster counts, so a redundant pick still recovers it
#' redundant <- which(sim$manifest$role == "redundant")[1]
#' score_selection_clusterwise(redundant, sim$manifest)$sensitivity
#' @export
score_selection_clusterwise <- function(selected, manifest) {
  stopifnot(
    "`manifest` must have `role` and `cluster_id` columns" =
      is.data.frame(manifest) && all(c("role", "cluster_id") %in% names(manifest))
  )
  p <- nrow(manifest)

  if (is.logical(selected)) {
    stopifnot("logical `selected` must have length nrow(manifest)" = length(selected) == p)
    sel_idx <- which(selected)
  } else {
    stopifnot(
      "`selected` must be feature indices or a logical vector" =
        is.numeric(selected) && all(is.finite(selected)) && all(selected == as.integer(selected)),
      "`selected` indices must lie in 1..nrow(manifest)" =
        length(selected) == 0L || (min(selected) >= 1L && max(selected) <= p)
    )
    sel_idx <- unique(as.integer(selected))
  }
  is_sel <- logical(p)
  is_sel[sel_idx] <- TRUE

  role <- manifest$role
  cid <- manifest$cluster_id
  rel_clusters <- unique(cid[role == "driver"])
  k <- length(rel_clusters)
  cluster_members <- function(cl) which(cid == cl & role %in% c("driver", "redundant"))

  n_recovered <- sum(vapply(rel_clusters, function(cl) any(is_sel[cluster_members(cl)]), logical(1)))
  sensitivity <- if (k > 0L) n_recovered / k else NA_real_

  fp <- sum(is_sel & role %in% c("distractor", "noise"))
  n_selected <- length(sel_idx)
  fdr <- if (n_selected > 0L) fp / n_selected else 0

  tp <- n_recovered
  recall <- sensitivity
  precision <- if (tp + fp > 0L) tp / (tp + fp) else NA_real_
  f1 <- if (isTRUE(precision > 0) && isTRUE(recall > 0)) {
    2 * precision * recall / (precision + recall)
  } else {
    0
  }

  exact_recovery <- (n_recovered == k) && (fp == 0L)

  parsimony_penalty <- sum(vapply(rel_clusters, function(cl) {
    max(0L, sum(is_sel[cluster_members(cl)]) - 1L)
  }, integer(1)))

  representative_per_cluster <- list()
  for (cl in rel_clusters) {
    members <- cluster_members(cl)
    sel_members <- members[is_sel[members]]
    if (length(sel_members) > 0L) {
      driver_member <- members[role[members] == "driver"]
      representative_per_cluster[[as.character(cl)]] <-
        c(intersect(sel_members, driver_member), setdiff(sel_members, driver_member))
    }
  }

  list(
    sensitivity = sensitivity,
    fdr = fdr,
    f1 = f1,
    exact_recovery = exact_recovery,
    parsimony_penalty = parsimony_penalty,
    representative_per_cluster = representative_per_cluster
  )
}
