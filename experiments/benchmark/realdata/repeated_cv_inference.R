# Inference over the repeated cross-validation of repeated_cv.R. Each cohort's score is the
# mean over its ten five-fold splits, which removes most of the resampling noise of a single
# split. Reports, for top-K at q = 2 against each rival at K = 10, the paired mean
# difference over cohorts with its standard error, the cohorts favouring each side and the
# Bayesian signed-rank posterior (argmax summary, ROPE 0.02 on Uno's C, 0.05 on stability),
# for Uno's C untruncated and truncated and for stability. It also reports how often the
# same comparison is decided when each repetition is analysed on its own, and the share of
# the variance of a paired difference that is between splits rather than between cohorts.
bench <- Sys.getenv("BENCH_DIR", ""); if (nzchar(bench) && dir.exists(bench)) setwd(bench)
if (requireNamespace("survfcbf", quietly = TRUE) &&
    "bayes_signed_rank" %in% getNamespaceExports("survfcbf")) {
  suppressMessages(library(survfcbf))
} else {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))
}
set.seed(1)
## RCV_IN / RCV_FOCUS select the run: the matched-size panels of repeated_cv.R (default) or
## the headline comparison of repeated_cv_main.R (RCV_IN=results/realdata_repeated_cv_main,
## RCV_FOCUS=survfcbf_topk)
IN <- Sys.getenv("RCV_IN", "results/realdata_repeated_cv")
OUTFILE <- file.path(IN, "inference.csv")
d <- do.call(rbind, lapply(list.files(IN, pattern = "^[A-Z]+_rep[0-9]+[.]csv$", full.names = TRUE),
                           utils::read.csv, check.names = FALSE))
stopifnot("expected 11 cohorts x 10 repetitions" = nrow(unique(d[, c("cohort", "repetition")])) == 110L,
          "repeated CV lost folds" = all(d$n_failed == 0))
FOCUS <- Sys.getenv("RCV_FOCUS", "survfcbf_q2_K10")
ROPES <- c(uno_c = 0.02, uno_c_trunc = 0.02, stability = 0.05)
verdict <- function(b) if (b[["p_x"]] > 0.95) "survfcbf better" else if (b[["p_y"]] > 0.95) "rival better" else
  if (b[["p_rope"]] > 0.95) "practically equivalent" else "undecided at n=11"
rows <- list()
for (rv in setdiff(unique(d$selector), FOCUS)) for (m in names(ROPES)) {
  w <- reshape(d[d$selector %in% c(FOCUS, rv), c("cohort", "repetition", "selector", m)],
               idvar = c("cohort", "repetition"), timevar = "selector", direction = "wide")
  w$diff <- w[[paste0(m, ".", FOCUS)]] - w[[paste0(m, ".", rv)]]
  cm <- aggregate(cbind(x = w[[paste0(m, ".", FOCUS)]], y = w[[paste0(m, ".", rv)]], diff = w$diff),
                  list(cohort = w$cohort), mean)
  b <- bayes_signed_rank(cm$x, cm$y, rope = ROPES[[m]])
  per_rep <- vapply(sort(unique(w$repetition)), function(r) {
    z <- w[w$repetition == r, ]
    verdict(bayes_signed_rank(z[[paste0(m, ".", FOCUS)]], z[[paste0(m, ".", rv)]], rope = ROPES[[m]], nsim = 20000))
  }, "")
  ## variance of the paired difference: between cohorts (of the cohort means) and within
  ## cohorts between splits
  within <- mean(tapply(w$diff, w$cohort, stats::var))
  rows[[length(rows) + 1L]] <- data.frame(
    rival = rv, metric = m, n_cohorts = nrow(cm), focus_mean = mean(cm$x), rival_mean = mean(cm$y),
    mean_diff = mean(cm$diff), se_diff = stats::sd(cm$diff) / sqrt(nrow(cm)),
    wins_focus = sum(cm$diff > 0), wins_rival = sum(cm$diff < 0),
    p_focus = b[["p_x"]], p_rope = b[["p_rope"]], p_rival = b[["p_y"]], verdict = verdict(b),
    reps_decided_same = sum(per_rep == verdict(b)), reps = length(per_rep),
    sd_between_splits = sqrt(within), sd_between_cohorts = stats::sd(cm$diff))
}
R <- do.call(rbind, rows)
utils::write.csv(R, OUTFILE, row.names = FALSE)
out <- R; num <- vapply(out, is.numeric, TRUE); out[num] <- lapply(out[num], round, 3)
print(out[, c("rival", "metric", "focus_mean", "rival_mean", "mean_diff", "se_diff", "wins_focus",
              "wins_rival", "p_focus", "p_rope", "p_rival", "verdict", "reps_decided_same",
              "sd_between_splits", "sd_between_cohorts")], row.names = FALSE)

## Sensitivity of the cohort-mean verdicts to the width of the ROPE and to the summary of
## the posterior (argmax, the reference convention used above, or the mean mass of each
## region), written apart so that inference.csv above does not change (it runs after the
## main loop, on its own draws)
RS <- list(uno_c_trunc = c(0.01, 0.02, 0.03), stability = c(0.03, 0.05, 0.07))
rs <- list()
for (rv in setdiff(unique(d$selector), FOCUS)) for (m in names(RS)) {
  w <- reshape(d[d$selector %in% c(FOCUS, rv), c("cohort", "repetition", "selector", m)],
               idvar = c("cohort", "repetition"), timevar = "selector", direction = "wide")
  cm <- aggregate(cbind(x = w[[paste0(m, ".", FOCUS)]], y = w[[paste0(m, ".", rv)]]), list(cohort = w$cohort), mean)
  for (rope in RS[[m]]) for (summ in c("argmax", "mean")) {
    b <- bayes_signed_rank(cm$x, cm$y, rope = rope, summary = summ)
    rs[[length(rs) + 1L]] <- data.frame(rival = rv, metric = m, rope = rope, summary = summ, p_focus = b[["p_x"]],
                                        p_rope = b[["p_rope"]], p_rival = b[["p_y"]], verdict = verdict(b))
  }
}
utils::write.csv(do.call(rbind, rs), file.path(IN, "rope_sensitivity.csv"), row.names = FALSE)
cat("REPEATED_CV_INFERENCE_DONE\n")
