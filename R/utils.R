# Internal helpers (not exported).

# Evaluate `expr` under a temporary RNG seed, restoring the global RNG state on
# exit. With `seed = NULL` the expression is evaluated against the current RNG
# state, so a caller that has set up an L'Ecuyer-CMRG stream (the orchestration
# layer) is respected and not overridden.
with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
  } else {
    on.exit(
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    )
  }
  set.seed(seed)
  expr
}

# Evaluate `expr` with `.Random.seed` set to a full RNG state vector (e.g. one
# L'Ecuyer-CMRG stream from parallel::nextRNGStream), restoring the previous
# global state on exit. Assigning a stream state switches the RNG kind
# implicitly, and restoring the saved vector switches it back, so no explicit
# RNGkind bookkeeping is needed. `state = NULL` evaluates against the current
# state, mirroring with_seed().
with_rng_state <- function(state, expr) {
  if (is.null(state)) {
    return(expr)
  }
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
  } else {
    on.exit(
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    )
  }
  assign(".Random.seed", state, envir = .GlobalEnv)
  expr
}

# Independent per-fold RNG stream states: seed a L'Ecuyer-CMRG generator from
# `seed` and space `n` streams with parallel::nextRNGStream(), restoring the
# caller's RNG kind and state afterwards. This is the package-wide seeding
# invariant (the simulator's run_replicates uses the same construction).
fold_rng_streams <- function(seed, n) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  old_kind <- RNGkind()
  ## Order matters and must stay in one handler: RNGkind() re-initialises the
  ## generator, so restoring the kind *after* the seed would overwrite the seed.
  on.exit({
    do.call(RNGkind, as.list(old_kind[1:2]))
    if (had_seed) {
      assign(".Random.seed", old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(as.integer(seed), kind = "L'Ecuyer-CMRG")
  st <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  streams <- vector("list", n)
  for (k in seq_len(n)) {
    streams[[k]] <- st
    st <- parallel::nextRNGStream(st)
  }
  streams
}

# Map a matrix of copula probabilities `U` (n x p, in (0, 1)) through per-column
# empirical inverse cdfs given `qcols` (a samples x p matrix whose column j holds
# the *sorted* reference values for feature j). Uses type-7 linear interpolation
# on the sorted reference, vectorised within each column.
empirical_quantile_matrix <- function(U, qcols) {
  m <- nrow(qcols)
  p <- ncol(qcols)
  n <- nrow(U)
  X <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    r <- qcols[, j]
    h <- (m - 1) * U[, j] + 1
    lo <- pmin(pmax(floor(h), 1), m)
    hi <- pmin(lo + 1, m)
    frac <- h - lo
    X[, j] <- r[lo] + frac * (r[hi] - r[lo])
  }
  X
}
