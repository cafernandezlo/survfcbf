# Cross-dataset comparison of the selectors over the 11 Bommert cohorts.
# Two complementary analyses, both implemented
# here with no scmamp/PMCMRplus dependency so the pipeline stays self-contained:
#
#   (1) Frequentist: average ranks + Friedman (Iman-Davenport F) + Nemenyi critical
#       difference, drawn as Demsar critical-difference diagrams.
#   (2) Bayesian: the Benavoli et al. (2017) Bayesian signed-rank test with a ROPE,
#       for the decision-relevant pairwise claims (is survfcbf_topk predictively
#       within a practical-equivalence margin of the baselines? how much worse on
#       stability?).
#
# A selector that returns NA on a cohort FAILED to produce a fittable model there;
# for the prediction ranks that counts as the worst outcome (failure is penalised,
# not dropped). Stability and parsimony are populated even on failed folds, so their
# ranks use all 11 cohorts. ROPE on the C-index is the pre-registered predictive-
# parity margin (0.02). Run in the container or locally with the package installed.

suppressWarnings(suppressMessages({
  bench <- Sys.getenv("BENCH_DIR", unset = ".")
  if (nzchar(bench) && dir.exists(bench)) setwd(bench)
}))
## bayes_signed_rank() lives in the package (one implementation for every analysis);
## fall back to the sources when the installed build predates it
if (requireNamespace("survfcbf", quietly = TRUE) &&
    "bayes_signed_rank" %in% getNamespaceExports("survfcbf")) {
  suppressMessages(library(survfcbf))
} else {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))
}
set.seed(1)

RESULTS <- "results/benchmark_realdata.csv"
OUTDIR  <- "results/realdata_crossdataset"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
ROPE_C    <- 0.02     # pre-registered predictive-parity margin on the C-index
ROPE_STAB <- 0.05     # practical-equivalence margin on the stability index
ALPHA     <- 0.05

d <- utils::read.csv(RESULTS, check.names = FALSE, stringsAsFactors = FALSE)
selectors <- unique(d$selector)
cohorts   <- unique(d$dataset)
FOCUS <- "survfcbf_topk"

## reshape one metric into a cohorts x selectors matrix ----------------------------
metric_matrix <- function(metric) {
  M <- matrix(NA_real_, length(cohorts), length(selectors),
              dimnames = list(cohorts, selectors))
  for (i in seq_len(nrow(d))) M[d$dataset[i], d$selector[i]] <- d[[metric]][i]
  M
}

## average ranks (rank 1 = best); NA ranked worst -------------------------------
avg_ranks <- function(M, higher_better = TRUE) {
  R <- t(apply(M, 1L, function(row) {
    v <- if (higher_better) row else -row
    v[is.na(v)] <- -Inf                                  # failure -> worst
    rank(-v, ties.method = "average")
  }))
  colMeans(R)
}

## Friedman (Iman-Davenport F) + Nemenyi critical difference --------------------
friedman_nemenyi <- function(M, higher_better = TRUE, alpha = ALPHA) {
  N <- nrow(M); k <- ncol(M)
  r <- avg_ranks(M, higher_better)
  chi2 <- 12 * N / (k * (k + 1)) * (sum(r^2) - k * (k + 1)^2 / 4)
  Ff   <- (N - 1) * chi2 / (N * (k - 1) - chi2)
  pval <- stats::pf(Ff, k - 1, (k - 1) * (N - 1), lower.tail = FALSE)
  q    <- stats::qtukey(1 - alpha, k, Inf) / sqrt(2)     # Nemenyi q_alpha
  CD   <- q * sqrt(k * (k + 1) / (6 * N))
  list(ranks = r, chi2 = chi2, F = Ff, p = pval, CD = CD, N = N, k = k)
}

## Demsar critical-difference diagram ------------------------------------------
plot_cd <- function(fn, ranks, CD, title) {
  k <- length(ranks); r <- sort(ranks)                   # rank 1 (best) on the left
  lo <- floor(min(r)); hi <- ceiling(max(r))
  grDevices::png(fn, width = 1100, height = 520, res = 130)
  on.exit(grDevices::dev.off())
  op <- graphics::par(mar = c(1, 1, 3, 1)); on.exit(graphics::par(op), add = TRUE)
  graphics::plot.new(); graphics::plot.window(xlim = c(lo - 0.5, hi + 0.5), ylim = c(-k - 1, 3))
  graphics::title(main = title)
  graphics::axis(3, at = lo:hi, pos = 0)
  ## stems: best ranks to the left, split half up / half down for label room
  ord <- order(r); half <- ceiling(k / 2)
  for (idx in seq_len(k)) {
    nm <- names(r)[idx]; xr <- r[idx]
    side_left <- idx <= half
    yend <- if (side_left) -idx else -(k - idx + 1)
    lab_x <- if (side_left) lo - 0.5 else hi + 0.5
    graphics::segments(xr, 0, xr, yend, col = "grey40")
    graphics::segments(xr, yend, lab_x, yend, col = "grey40")
    graphics::text(lab_x, yend, labels = sprintf("%s  (%.2f)", nm, xr),
                   pos = if (side_left) 4 else 2, xpd = TRUE,
                   font = if (nm == FOCUS) 2 else 1,
                   col  = if (nm == FOCUS) "firebrick" else "black")
  }
  ## CD bar
  graphics::segments(lo, 2, lo + CD, 2, lwd = 3)
  graphics::text(lo + CD / 2, 2.4, sprintf("CD = %.2f", CD), cex = 0.9)
  ## cliques: connect groups whose rank spread <= CD
  rs <- sort(r); y <- -0.4; done_hi <- -Inf
  for (i in seq_len(k)) {
    j <- max(which(rs <= rs[i] + CD))
    if (j > i && rs[j] > done_hi) {
      graphics::segments(rs[i] - 0.03, y, rs[j] + 0.03, y, lwd = 4, col = "steelblue")
      y <- y - 0.35; done_hi <- rs[j]
    }
  }
  invisible(fn)
}

# ============================ run ============================================
metrics <- list(
  list(key = "uno_c",     hi = TRUE,  lab = "Uno's C (discrimination)"),
  list(key = "ibs",       hi = FALSE, lab = "Integrated Brier score"),
  list(key = "stability", hi = TRUE,  lab = "Stability (SMA-Count)"),
  list(key = "parsimony", hi = FALSE, lab = "Parsimony (n features)")
)

cat("=== validity (cohorts with a fittable model / mean failed folds) ===\n")
valid <- data.frame(selector = selectors,
  na_uno   = sapply(selectors, function(s) sum(is.na(d$uno_c[d$selector == s]))),
  mean_failed = round(sapply(selectors, function(s) mean(d$n_failed[d$selector == s])), 2),
  row.names = NULL)
print(valid[order(valid$na_uno), ], row.names = FALSE)

cat("\n=== Friedman + Nemenyi (avg ranks, rank 1 = best) ===\n")
summ <- list()
for (mt in metrics) {
  M  <- metric_matrix(mt$key)
  fn <- friedman_nemenyi(M, higher_better = mt$hi)
  summ[[mt$key]] <- fn
  cat(sprintf("\n-- %s --  Friedman F=%.2f, p=%.3g, CD=%.2f (N=%d, k=%d)\n",
              mt$lab, fn$F, fn$p, fn$CD, fn$N, fn$k))
  rr <- sort(fn$ranks)
  for (nm in names(rr)) cat(sprintf("   %-16s %.2f%s\n", nm, rr[nm],
                                    if (nm == FOCUS) "   <-- survfcbf_topk" else ""))
  plot_cd(file.path(OUTDIR, sprintf("cd_%s.png", mt$key)), fn$ranks, fn$CD, mt$lab)
}

cat("\n=== Bayesian signed-rank: survfcbf_topk vs each rival ===\n")
rivals <- setdiff(selectors, FOCUS)
Mc <- metric_matrix("uno_c"); Ms <- metric_matrix("stability")
bt <- data.frame()
## summary = "argmax" (Benavoli reference convention). The package's mean-mass default is
## bounded below 1 at this effect-to-heterogeneity ratio and hid one decided pair in each
## direction; under argmax a region exceeding 0.95 is a decision in the standard sense.
## The mean-mass summary is reported alongside (C_*_mean), because it decides far fewer
## pairs and a reader should see both.
for (rv in rivals) {
  bc <- bayes_signed_rank(Mc[, FOCUS], Mc[, rv], rope = ROPE_C,   higher_better = TRUE, summary = "argmax")
  bs <- bayes_signed_rank(Ms[, FOCUS], Ms[, rv], rope = ROPE_STAB, higher_better = TRUE, summary = "argmax")
  bm <- bayes_signed_rank(Mc[, FOCUS], Mc[, rv], rope = ROPE_C,   higher_better = TRUE, summary = "mean")
  bt <- rbind(bt, data.frame(rival = rv,
    C_topk = round(bc["p_x"], 2), C_rope = round(bc["p_rope"], 2), C_rival = round(bc["p_y"], 2), C_n = bc["n"],
    S_topk = round(bs["p_x"], 2), S_rope = round(bs["p_rope"], 2), S_rival = round(bs["p_y"], 2),
    S_n = bs["n"],
    C_topk_mean = round(bm["p_x"], 2), C_rope_mean = round(bm["p_rope"], 2), C_rival_mean = round(bm["p_y"], 2)))
}
rownames(bt) <- NULL
cat(sprintf("(C-index ROPE +/-%.2f, stability ROPE +/-%.2f; probs = P(topk better)/P(rope)/P(rival better))\n",
            ROPE_C, ROPE_STAB))
print(bt, row.names = FALSE)

saveRDS(list(valid = valid, friedman = summ, bayes = bt), file.path(OUTDIR, "crossdataset.rds"))
utils::write.csv(bt, file.path(OUTDIR, "bayes_signed_rank.csv"), row.names = FALSE)
cat(sprintf("\nfigures + tables -> %s\n", OUTDIR))
cat("CROSSDATASET_DONE\n")
