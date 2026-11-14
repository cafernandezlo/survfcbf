#' Run Monte Carlo replicates of a scenario with per-replicate L'Ecuyer streams
#'
#' Repeatedly evaluates [simulate_dataset()] for a fixed scenario, seeding each
#' replicate with its own L'Ecuyer-CMRG stream derived from a single
#' `master_seed`, and summarises the realised diagnostics with their Monte Carlo
#' standard errors. A sweep is this function applied to a set
#' of one-factor-at-a-time scenarios.
#'
#' @details
#' Independent, reproducible streams are spaced with [parallel::nextRNGStream()]
#' so the replicates do not share random numbers and the whole run reproduces
#' from `master_seed`. Each replicate is run with the scenario's `seed` forced to
#' `NULL` so it draws from the stream set here. The global RNG kind and state are
#' restored on exit.
#'
#' @param params A scenario config list, normally from [default_params()].
#' @param R Integer. Number of Monte Carlo replicates. Default `200`.
#' @param master_seed Integer master seed from which the per-replicate streams
#'   are derived. Default `1`.
#'
#' @return A list with `summary` (a data.frame of `diagnostic`, `mean` and
#'   `mc_se = sd / sqrt(R)` for each realised diagnostic), `replicates` (the
#'   diagnostics-by-replicate matrix), `R`, `master_seed`, and `params`.
#'
#' @seealso [simulate_dataset()], [default_params()]
#' @export
run_replicates <- function(params, R = 200L, master_seed = 1L) {
  stopifnot(
    "`params` must be a list"                  = is.list(params),
    "`R` must be a single positive integer"    = length(R) == 1L && is.numeric(R) && R >= 1,
    "`master_seed` must be a single number"    = length(master_seed) == 1L && is.numeric(master_seed)
  )
  R <- as.integer(R)

  old_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })

  RNGkind("L'Ecuyer-CMRG")
  set.seed(master_seed)
  stream <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  pp <- params
  pp$seed <- NULL                       # streams are managed here, not inside simulate_dataset

  keys <- c("beta_scale", "censoring_fraction", "n_events", "concordance", "intracluster_cor")
  mat <- matrix(NA_real_, nrow = length(keys), ncol = R, dimnames = list(keys, NULL))
  for (r in seq_len(R)) {
    assign(".Random.seed", stream, envir = .GlobalEnv)
    realized <- simulate_dataset(pp)$realized
    mat[, r] <- unlist(realized[keys])
    stream <- parallel::nextRNGStream(stream)
  }

  summary <- data.frame(
    diagnostic = keys,
    mean = rowMeans(mat),
    mc_se = apply(mat, 1L, stats::sd) / sqrt(R),
    row.names = NULL
  )
  list(summary = summary, replicates = mat, R = R, master_seed = master_seed, params = params)
}
