# Does doubling the outer folds shrink the paired between-cohort spread?
# An empirical test of heterogeneity versus estimation noise.
#
# If the spread of paired Uno's C differences were pure between-cohort heterogeneity,
# doubling the folds would leave it unchanged (ratio ~1). If it were pure evaluation
# noise, doubling the folds would shrink it towards 1/sqrt(2) ~ 0.71. The observed ratio
# lands in between and is itself estimated from only 11 paired differences per selector
# pair, so it carries large sampling error: per-pair ratios range well past both anchors.
# The honest conclusion is therefore directional, not quantitative: evaluation noise is a
# material component of the spread, more folds do buy back some resolution, and the
# required-datasets figures remain valid FOR THE DESIGN THEY DESCRIBE (they use the
# variability that a study with that protocol actually faces).
#
# Usage (after the 10-fold array, folds_sensitivity.R, has produced results/realdata_folds):
#   Rscript --no-init-file realdata/folds_compare.R
bench <- Sys.getenv("BENCH_DIR", ""); if (nzchar(bench) && dir.exists(bench)) setwd(bench)
OUT <- "results/realdata_folds"

f10 <- do.call(rbind, lapply(list.files(OUT, pattern = "_f10\\.csv$", full.names = TRUE),
                             utils::read.csv, check.names = FALSE))
cf <- list.files("results/realdata_cardinality", pattern = "^[A-Z]+\\.csv$", full.names = TRUE)
f5 <- do.call(rbind, lapply(cf, utils::read.csv, check.names = FALSE))
f5 <- f5[f5$K == 20, ]

pair_sd <- function(d, selcol) {
  sels <- sort(unique(d[[selcol]])); coh <- sort(unique(d$cohort)); out <- c()
  for (i in seq_along(sels)) for (j in seq_along(sels)) if (i < j) {
    x <- setNames(rep(NA_real_, length(coh)), coh); y <- x
    zx <- d[d[[selcol]] == sels[i], ]; zy <- d[d[[selcol]] == sels[j], ]
    x[zx$cohort] <- zx$uno_c; y[zy$cohort] <- zy$uno_c
    dd <- (x - y)[is.finite(x) & is.finite(y)]
    if (length(dd) >= 5) out <- c(out, stats::setNames(stats::sd(dd), paste(sels[i], "vs", sels[j])))
  }
  out
}
s5 <- pair_sd(f5, "base"); s10 <- pair_sd(f10, "selector")
common <- intersect(names(s5), names(s10))
res <- data.frame(pair = common, sd_5fold = unname(s5[common]), sd_10fold = unname(s10[common]),
                  ratio = unname(s10[common] / s5[common]))
utils::write.csv(res, file.path(OUT, "folds_comparison.csv"), row.names = FALSE)

cat("=== paired SD of Uno's C differences at K=20, 5 versus 10 outer folds ===\n")
print(within(res, { sd_5fold <- round(sd_5fold, 4); sd_10fold <- round(sd_10fold, 4)
                    ratio <- round(ratio, 2) }), row.names = FALSE)
cat(sprintf("\nmedians: 5-fold %.4f | 10-fold %.4f | ratio %.2f\n",
            stats::median(res$sd_5fold), stats::median(res$sd_10fold), stats::median(res$ratio)))
cat("anchors: 1.00 = pure heterogeneity, 0.71 = pure evaluation noise\n")
n_needed <- function(delta, sigma) ceiling((stats::qnorm(0.975) + stats::qnorm(0.8))^2 * sigma^2 / delta^2)
cat(sprintf("datasets for delta=0.03 at this K: %d with 5 folds, %d with 10 folds\n",
            n_needed(0.03, stats::median(res$sd_5fold)), n_needed(0.03, stats::median(res$sd_10fold))))
cat("FOLDS_COMPARE_DONE\n")
