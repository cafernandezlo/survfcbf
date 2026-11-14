# Cross-validation of the information-theoretic core against the reference FCBF package
# (the correctness claims of the Application Note). The checks were first run by hand and
# are consolidated here with the same seeds and the same order of random draws, so the
# counts reproduce rather than approximate. They are therefore plain set.seed()
# sweeps and not L'Ecuyer streams: each one is a fixed list of test cases, not a Monte Carlo
# estimate. Output: results/validation/fcbf_crossvalidation.csv, read by make_statistics.R.
#
# Needs FCBF, which left Bioconductor with release 3.20:
#   remotes::install_github("lubianat/FCBF")   # R/ identical to the last Bioconductor release
# and SummarizedExperiment + SingleCellExperiment for its scDengue dataset.
#
# Usage (from experiments/benchmark):  Rscript --no-init-file realdata/fcbf_crossvalidation.R
setwd(Sys.getenv("BENCH_DIR", unset = "."))
suppressMessages(devtools::load_all("../..", quiet = TRUE))
for (pkg in c("FCBF", "SummarizedExperiment", "SingleCellExperiment"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("this script needs the package ", pkg, call. = FALSE)
suppressMessages(library(FCBF))
TIE_TOL <- 1e-12

row <- function(check, comparisons, mismatches, skipped_ties = 0L, mismatches_at_tie = 0L,
                max_abs_diff = NA_real_, detail = "")
  data.frame(check, comparisons, mismatches, skipped_ties, mismatches_at_tie, max_abs_diff, detail)
report <- function(r) { cat(sprintf("%-26s comparisons %5d | mismatches %3d | skipped ties %2d | at a tie %3d | max diff %s | %s\n",
                                    r$check, r$comparisons, r$mismatches, r$skipped_ties, r$mismatches_at_tie,
                                    format(r$max_abs_diff, digits = 3), r$detail)); r }

## the reference, silenced, as integer column indices in relevance order (NULL on error)
reference <- function(X, y, delta, samples_in_rows = TRUE) {
  out <- NULL
  ok <- tryCatch({ utils::capture.output(
      out <- FCBF::fcbf(X, as.factor(y), minimum_su = delta, samples_in_rows = samples_in_rows,
                        verbose = FALSE)); TRUE },
    error = function(e) FALSE)
  if (ok) as.integer(out$index) else NULL
}
## TRUE when some SU(Fi, Fj) equals SU(Fj, y): the one place where >= and > can disagree
has_exact_tie <- function(X, y) {
  p <- ncol(X)
  suy <- vapply(X, function(z) symmetrical_uncertainty(z, y), numeric(1))
  M <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) if (i != j) M[i, j] <- symmetrical_uncertainty(X[[i]], X[[j]])
  any(vapply(seq_len(p), function(j) any(abs(M[-j, j] - suy[j]) < TIE_TOL), logical(1)))
}

## ---------------------------------------------------------------- symmetrical uncertainty
su_pairs <- function() {
  set.seed(1)
  worst <- 0; n_pairs <- 0L
  for (rep in 1:200) {
    n  <- sample(20:120, 1); kx <- sample(2:6, 1); ky <- sample(2:6, 1)
    x <- sample.int(kx, n, replace = TRUE); y <- sample.int(ky, n, replace = TRUE)
    b <- FCBF::get_SU_for_vector_pair(as.factor(x), as.factor(y))
    if (is.nan(b)) next
    worst <- max(worst, abs(symmetrical_uncertainty(x, y) - b)); n_pairs <- n_pairs + 1L
  }
  row("su_random_pairs", n_pairs, 0L, max_abs_diff = worst, detail = "n 20-120, 2-6 levels")
}

## ---------------------------------------------------------------- selections, no exact ties
selections_nondegenerate <- function() {
  set.seed(2024)
  n_run <- 0L; n_mis <- 0L; n_tie <- 0L
  for (rep in 1:400) {
    n  <- sample(60:250, 1); p <- sample(5:14, 1); ky <- sample(2:4, 1)
    y  <- sample.int(ky, n, replace = TRUE)
    kx <- sample(2:5, 1)
    cols <- lapply(seq_len(p), function(j) {
      pr <- runif(1, 0.15, 0.85)                            # always noisy: never equals y
      ifelse(runif(n) < pr, sample.int(kx, n, TRUE), ((y - 1) %% kx) + 1)
    })
    for (k in seq_len(sample(0:3, 1))) {                    # redundant clusters: noisy copies
      src <- sample(seq_len(p), 1)
      cols[[sample(seq_len(p), 1)]] <-
        ifelse(runif(n) < runif(1, 0.05, 0.4), sample.int(kx, n, TRUE), cols[[src]])
    }
    X <- as.data.frame(lapply(cols, function(z) factor(letters[z])))
    names(X) <- paste0("V", seq_len(p))
    keep <- vapply(X, function(z) length(unique(z)) > 1, logical(1))  # constant: undefined SU there
    X <- X[, keep, drop = FALSE]
    if (ncol(X) < 3) next
    if (has_exact_tie(X, y)) { n_tie <- n_tie + 1L; next }
    delta <- sample(c(0, 0.01, 0.03, 0.06), 1)
    a <- as.integer(fcbf_select(X, y, delta = delta)$selected)
    b <- reference(X, y, delta)
    if (is.null(b)) next
    n_run <- n_run + 1L
    if (!identical(a, b)) n_mis <- n_mis + 1L
  }
  row("selections_nondegenerate", n_run, n_mis, skipped_ties = n_tie, detail = "400 generated, same features in the same order")
}

## ---------------------------------------------------------------- selections, ties allowed
## Exact duplicates and copies of the target are allowed here, so the two rules can disagree.
## Every disagreement has to sit at an exact tie, which is the documented divergence.
selections_with_ties <- function() {
  set.seed(7)
  n_run <- 0L; n_mis <- 0L; n_at_tie <- 0L
  for (rep in 1:300) {
    n  <- sample(40:200, 1); p <- sample(4:12, 1); ky <- sample(2:4, 1)
    y  <- sample.int(ky, n, replace = TRUE)
    kx <- sample(2:4, 1)
    base <- lapply(seq_len(p), function(j) {
      pr <- runif(1, 0, 0.9)
      ifelse(runif(n) < pr, sample.int(kx, n, TRUE), ((y - 1) %% kx) + 1)
    })
    if (runif(1) < 0.5) {
      src <- sample(seq_len(p), 1); dst <- sample(setdiff(seq_len(p), src), 1)
      base[[dst]] <- base[[src]]
    }
    X <- as.data.frame(lapply(base, function(z) factor(letters[z])))
    names(X) <- paste0("V", seq_len(p))
    delta <- sample(c(0, 0.01, 0.03, 0.08), 1)
    a <- as.integer(fcbf_select(X, y, delta = delta)$selected)
    b <- reference(X, y, delta)
    if (is.null(b)) next
    n_run <- n_run + 1L
    if (!identical(a, b)) {
      n_mis <- n_mis + 1L
      keep <- vapply(X, function(z) length(unique(z)) > 1, logical(1))
      if (has_exact_tie(X[, keep, drop = FALSE], y)) n_at_tie <- n_at_tie + 1L
    }
  }
  row("selections_with_ties", n_run, n_mis, mismatches_at_tie = n_at_tie, detail = "300 generated, duplicates allowed")
}

## ---------------------------------------------------------------- top-K is an exact early stop
top_vs_head <- function(X, y, delta, Ks, ...) {
  full <- as.integer(fcbf_select(X, y, delta = delta, ...)$selected)
  vapply(Ks, function(K) !identical(as.integer(fcbf_select(X, y, delta = delta, top = K, ...)$selected),
                                    utils::head(full, K)), logical(1))
}
topk_random <- function() {
  set.seed(99)
  bad <- logical()
  for (rep in 1:400) {
    n  <- sample(40:200, 1); p <- sample(4:20, 1); ky <- sample(2:5, 1)
    y  <- sample.int(ky, n, replace = TRUE)
    kx <- sample(2:5, 1)
    cols <- lapply(seq_len(p), function(j) {
      pr <- runif(1, 0, 0.95)
      ifelse(runif(n) < pr, sample.int(kx, n, TRUE), ((y - 1) %% kx) + 1)
    })
    for (k in seq_len(sample(0:4, 1))) {
      src <- sample(seq_len(p), 1)
      cols[[sample(seq_len(p), 1)]] <-
        ifelse(runif(n) < runif(1, 0, 0.5), sample.int(kx, n, TRUE), cols[[src]])
    }
    X <- as.data.frame(cols); names(X) <- paste0("V", seq_len(p))
    delta <- sample(c(0, 0.01, 0.05, 0.15), 1)
    bad <- c(bad, top_vs_head(X, y, delta, seq_len(min(p, 6))))
  }
  row("topk_random", length(bad), sum(bad), detail = "400 datasets, K 1-6, four deltas")
}
topk_relevance_weights <- function() {
  set.seed(5)
  bad <- logical()
  for (rep in 1:150) {
    n <- sample(50:150, 1); p <- sample(5:15, 1)
    y <- sample.int(3, n, TRUE)
    X <- as.data.frame(replicate(p, sample.int(3, n, TRUE)))
    if (p >= 3) X[[2]] <- X[[1]]                            # exact duplicate
    rel <- runif(p)
    w   <- runif(n, 0.5, 2)
    for (cfg in list(list(rel = rel, w = NULL), list(rel = NULL, w = w), list(rel = rel, w = w)))
      bad <- c(bad, top_vs_head(X, y, 0, 1:4, relevance = cfg$rel, weights = cfg$w))
  }
  row("topk_relevance_weights", length(bad), sum(bad), detail = "precomputed relevance and observation weights")
}
topk_ties <- function() {
  set.seed(777)
  bad <- logical()
  for (rep in 1:300) {
    n <- sample(30:120, 1); p <- sample(6:16, 1)
    y <- sample.int(sample(2:4, 1), n, TRUE)
    X <- as.data.frame(replicate(p, sample.int(sample(2:4, 1), n, TRUE)))
    for (k in 1:sample(1:5, 1)) X[[sample(p, 1)]] <- X[[sample(p, 1)]]    # many duplicated columns
    rel <- sample(c(0.9, 0.5, 0.1), p, TRUE)                # only 3 distinct relevance values
    for (cfg in list(NULL, rel))
      bad <- c(bad, top_vs_head(X, y, 0, seq_len(min(p, 5)), relevance = cfg))
  }
  row("topk_ties", length(bad), sum(bad), detail = "duplicated columns and tied relevance")
}

## ---------------------------------------------------------------- the reference's own data
## Discretised by the reference's own discretize_exprs, so discretisation cannot be a confounder.
scdengue <- function() {
  e <- new.env(); utils::data("scDengue", package = "FCBF", envir = e)
  exprs <- as.data.frame(SummarizedExperiment::assay(e$scDengue, "logcounts"))
  y <- as.factor(SummarizedExperiment::colData(e$scDengue)$infection)
  disc <- FCBF::discretize_exprs(exprs)                     # genes in rows, cells in columns
  Xs <- as.data.frame(t(disc), stringsAsFactors = FALSE); colnames(Xs) <- rownames(disc)
  deltas <- c(0.03, 0.05, 0.08, 0.1, 0.13)
  mine <- lapply(deltas, function(d) as.integer(fcbf_select(Xs, y, delta = d)$selected))
  theirs <- lapply(deltas, function(d) reference(disc, y, d, samples_in_rows = FALSE))
  same <- mapply(identical, mine, theirs)
  bad <- top_vs_head(Xs, y, 0.05, 1:10)
  list(row("selections_scdengue", length(deltas), sum(!same),
           detail = sprintf("%d genes x %d cells, panels of %s", nrow(disc), ncol(disc),
                            paste(lengths(mine), collapse = "/"))),
       row("topk_scdengue", length(bad), sum(bad), detail = "delta 0.05, K 1-10"))
}

res <- do.call(rbind, lapply(c(list(su_pairs(), selections_nondegenerate(), selections_with_ties(),
                                    topk_random(), topk_relevance_weights(), topk_ties()), scdengue()), report))
OUT <- "results/validation"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(OUT, "fcbf_crossvalidation.csv"), row.names = FALSE)
cat("\nsessionInfo FCBF:", as.character(utils::packageVersion("FCBF")), "\nFCBF_XVAL_DONE\n")
