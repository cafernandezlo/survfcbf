# Uncertainty for the cardinality-curve claims. The curve itself
# (cardinality_curve.R) reports per-cohort Uno's C; this attaches the inference, because
# "survfcbf beats variance" and "mim beats survfcbf" are the comparisons that matter, and a
# mean over 11 cohorts is not evidence on its own. For every K and every rival,
# reports the paired mean difference with its standard error, how many cohorts favour each
# side, and the Bayesian signed-rank posterior against the pre-registered 0.02 ROPE on the
# C-index ([bayes_signed_rank()], the same machinery as the cross-dataset analysis).
bench <- Sys.getenv("BENCH_DIR", ""); if (nzchar(bench) && dir.exists(bench)) setwd(bench)
## the installed package may predate bayes_signed_rank(); fall back to the sources
if (requireNamespace("survfcbf", quietly = TRUE) &&
    "bayes_signed_rank" %in% getNamespaceExports("survfcbf")) {
  suppressMessages(library(survfcbf))
} else {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))
}
set.seed(1)

ROPE <- 0.02
FOCUS <- "survfcbf"
IN  <- "results/realdata_cardinality"
OUT <- "results/realdata_cardinality"

## per-cohort cells only: this script writes its own CSV into the same directory, so it
## must not read it back on a re-run (that is what makes the stage idempotent)
OUTFILE <- "cardinality_inference.csv"
files <- setdiff(list.files(IN, pattern = "csv$", full.names = TRUE),
                 file.path(IN, OUTFILE))
stopifnot("no cardinality cells found" = length(files) > 0)
d <- do.call(rbind, lapply(files, utils::read.csv, check.names = FALSE))
cohorts <- sort(unique(d$cohort))
cat(sprintf("cohorts=%d  rows=%d  selectors=%s\n", length(cohorts), nrow(d),
            paste(sort(unique(d$base)), collapse = ", ")))

## per-cohort Uno's C matrix for one selector at one K (NA where it did not run)
score_of <- function(base, K) {
  v <- setNames(rep(NA_real_, length(cohorts)), cohorts)
  s <- d[d$base == base & d$K == K, ]
  v[s$cohort] <- s$uno_c
  v
}

rivals <- setdiff(unique(d$base), FOCUS)
rows <- list()
for (K in sort(unique(d$K[d$base == FOCUS]))) {
  x <- score_of(FOCUS, K)
  for (rv in rivals) {
    y <- score_of(rv, K)
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2L) next
    dif <- (x - y)[ok]
    ## summary = "argmax": the Benavoli reference convention (fraction of draws in which
    ## each region holds the largest mass). The package default, the mean mass per region,
    ## is bounded by P((Z1+Z2)/2 > rope) and stalls near 0.60 at this effect size, and a
    ## previous version of this script mistook that bound for a property of Bayesian
    ## comparison and leaned on "p_rival < 0.05" as the decision -- a test whose success
    ## rate FALLS with more datasets (p_rival converges to 0.099 at delta = 0.03). Under
    ## argmax the conventional 0.95 rule behaves normally, so it is used directly.
    b <- bayes_signed_rank(x, y, rope = ROPE, summary = "argmax")
    sdev <- stats::sd(dif)
    rows[[length(rows) + 1L]] <- data.frame(
      K = K, rival = rv, n_cohorts = sum(ok),
      mean_diff = mean(dif), se_diff = sdev / sqrt(sum(ok)),
      wins_focus = sum(dif > 0), wins_rival = sum(dif < 0),
      p_focus = b[["p_x"]], p_rope = b[["p_rope"]], p_rival = b[["p_y"]],
      stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, rows)
## Standard rule on the argmax posterior: declare whichever region exceeds 0.95,
## otherwise the 11 cohorts do not decide the pair.
R$verdict <- ifelse(R$p_focus > 0.95, "survfcbf better",
             ifelse(R$p_rival > 0.95, "rival better",
             ifelse(R$p_rope  > 0.95, "practically equivalent", "undecided at n=11")))
utils::write.csv(R, file.path(OUT, OUTFILE), row.names = FALSE)

fmt <- function(z) data.frame(K = z$K, rival = z$rival, n = z$n_cohorts,
  diff = sprintf("%+.3f +/- %.3f", z$mean_diff, z$se_diff),
  cohorts = sprintf("%d-%d", z$wins_focus, z$wins_rival),
  P_survfcbf = round(z$p_focus, 2), P_rope = round(z$p_rope, 2), P_rival = round(z$p_rival, 2),
  verdict = z$verdict)
cat(sprintf("\n=== survfcbf vs each rival, paired over cohorts (ROPE +/-%.2f on Uno's C) ===\n", ROPE))
for (rv in rivals) {
  z <- R[R$rival == rv, ]
  if (!nrow(z)) next
  cat(sprintf("\n-- vs %s --\n", rv))
  print(fmt(z[order(z$K), ]), row.names = FALSE)
}
cat("\n(diff = mean paired Uno's C difference +/- SE; cohorts = wins for survfcbf - wins for rival)\n")
cat("CARDINALITY_INFERENCE_DONE\n")
