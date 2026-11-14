# Inference for the mRMR comparison of panel_redundancy.R: survfcbf top-K (q = 2) against
# mRMR at K = 5, 10 and 20, paired over the eleven TCGA cohorts. For Uno's C (ROPE 0.02)
# and stability (ROPE 0.05), the paired mean difference with its standard error, the
# cohorts favouring each side and the Bayesian signed-rank posterior under the argmax
# summary, with the 0.95 decision rule of cardinality_inference.R.
bench <- Sys.getenv("BENCH_DIR", ""); if (nzchar(bench) && dir.exists(bench)) setwd(bench)
## the installed package may predate bayes_signed_rank(); fall back to the sources
if (requireNamespace("survfcbf", quietly = TRUE) &&
    "bayes_signed_rank" %in% getNamespaceExports("survfcbf")) {
  suppressMessages(library(survfcbf))
} else {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))
}
set.seed(1)

IN <- "results/realdata_redundancy"
OUTFILE <- file.path(IN, "mrmr_inference.csv")
ROPES <- c(uno_c = 0.02, stability = 0.05)
d <- do.call(rbind, lapply(list.files(IN, pattern = "^[A-Z]+\\.csv$", full.names = TRUE),
                           utils::read.csv, check.names = FALSE))
cohorts <- sort(unique(d$cohort))
stopifnot("expected the eleven TCGA cohorts" = length(cohorts) == 11L)
score_of <- function(sel, metric) { z <- d[d$selector == sel, ]; z[[metric]][match(cohorts, z$cohort)] }

rows <- list()
for (K in c(5L, 10L, 20L)) for (metric in names(ROPES)) {
  x <- score_of(sprintf("survfcbf_q2_K%d", K), metric)
  y <- score_of(sprintf("mrmr_K%d", K), metric)
  dif <- x - y
  b <- bayes_signed_rank(x, y, rope = ROPES[[metric]], summary = "argmax")
  rows[[length(rows) + 1L]] <- data.frame(
    K = K, metric = metric, rival = "mrmr", n_cohorts = sum(is.finite(dif)),
    mean_diff = mean(dif), se_diff = stats::sd(dif) / sqrt(sum(is.finite(dif))),
    wins_focus = sum(dif > 0), wins_rival = sum(dif < 0),
    p_focus = b[["p_x"]], p_rope = b[["p_rope"]], p_rival = b[["p_y"]])
}
R <- do.call(rbind, rows)
R$verdict <- ifelse(R$p_focus > 0.95, "survfcbf better",
             ifelse(R$p_rival > 0.95, "rival better",
             ifelse(R$p_rope  > 0.95, "practically equivalent", "undecided at n=11")))
utils::write.csv(R, OUTFILE, row.names = FALSE)
print(within(R, { mean_diff <- round(mean_diff, 3); se_diff <- round(se_diff, 3)
                  p_focus <- round(p_focus, 2); p_rope <- round(p_rope, 2); p_rival <- round(p_rival, 2) }),
      row.names = FALSE)
cat("REDUNDANCY_INFERENCE_DONE\n")
