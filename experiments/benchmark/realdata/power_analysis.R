# How many datasets does a survival feature-selection benchmark actually need?
# Every comparison on our 11 cohorts came
# back inconclusive, which raises the question the field never asks: what between-dataset
# variability do these benchmarks face, and how many datasets would a claim of "selector A
# beats selector B by delta" require?
#
# Two layers:
#   (1) analytic, paired t-test: n = (z_{1-a/2} + z_{1-b})^2 * sigma^2 / delta^2, the
#       textbook reference;
#   (2) simulation of the test we actually use -- the Bayesian signed-rank with a ROPE --
#       because with a ROPE of r a true difference below r is practically equivalent BY
#       DEFINITION and no sample size can declare it a win. That floor is a property of the
#       decision rule, not a failure of the benchmark, and has to be reported as such.
#
# sigma is estimated EMPIRICALLY from the paired per-cohort differences we observed, over
# all selector pairs, in the primary metric: Uno's C truncated at the 80th percentile of the
# training follow-up (realdata/uno_truncated.R), averaged over the ten five-fold splits of
# each cohort (realdata/repeated_cv_main.R for the headline comparison, realdata/repeated_cv.R
# for the matched sizes K = 5, 10, 20). sigma is the median SD within one of two families of
# comparisons, each with its own power statement:
#   headline  top-K (ten genes) against the variance filter, CAR and MIM at the ten-percent
#             cut-off, the comparison the note makes (3 pairs);
#   matched   top-K at q = 2 against the same three rankers cut to the same K = 5, 10, 20
#             genes (9 pairs).
# A median over every selector pair of a run (the earlier rule) depends on which variants
# the run happens to contain, since near-identical variants (mRMR at q = 2 against top-K)
# give tiny SDs, so the other pairs are kept in pairwise_sd.csv (family "other") but enter
# no sigma. The single-split, untruncated SDs of the original benchmark and cardinality
# sweep stay there too (family "legacy"), for context only.
bench <- Sys.getenv("BENCH_DIR", ""); if (nzchar(bench) && dir.exists(bench)) setwd(bench)
if (requireNamespace("survfcbf", quietly = TRUE) &&
    "bayes_signed_rank" %in% getNamespaceExports("survfcbf")) {
  suppressMessages(library(survfcbf))
} else {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))
}
set.seed(1)
OUT <- "results/realdata_power"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ROPE <- 0.02

## ---- 1. empirical between-dataset SD of paired differences ---------------------
pair_sds <- function(M) {                       # M: cohorts x selectors matrix of Uno's C
  sel <- colnames(M); out <- list()
  for (i in seq_along(sel)) for (j in seq_along(sel)) if (i < j) {
    d <- M[, i] - M[, j]; d <- d[is.finite(d)]
    if (length(d) >= 5L) out[[length(out) + 1L]] <- data.frame(
      a = sel[i], b = sel[j], n = length(d), mean_diff = mean(d), sd_diff = stats::sd(d))
  }
  do.call(rbind, out)
}
sds <- list()

## primary metric: cohort means over ten splits of the truncated Uno's C
RCV <- c(main = "results/realdata_repeated_cv_main", K5 = "results/realdata_repeated_cv_K5",
         K10 = "results/realdata_repeated_cv", K20 = "results/realdata_repeated_cv_K20")
## POWER_RCV_SETS restricts the runs used (e.g. "main,K10" for a provisional estimate)
RCV <- RCV[strsplit(Sys.getenv("POWER_RCV_SETS", paste(names(RCV), collapse = ",")), ",", fixed = TRUE)[[1]]]
stopifnot("unknown POWER_RCV_SETS entry" = !anyNA(RCV))
for (nm in names(RCV)) {
  f <- list.files(RCV[[nm]], pattern = "^[A-Z]+_rep[0-9]+[.]csv$", full.names = TRUE)
  stopifnot("missing repeated CV run" = length(f) > 0L)
  r <- do.call(rbind, lapply(f, utils::read.csv, check.names = FALSE))
  stopifnot("expected 11 cohorts x 10 repetitions" = nrow(unique(r[, c("cohort", "repetition")])) == 110L)
  M <- tapply(r$uno_c_trunc, list(r$cohort, r$selector), mean)
  sds[[nm]] <- cbind(source = paste0("rcv_", nm), pair_sds(M))
}

## legacy: one split, untruncated (context only)
b <- utils::read.csv("results/benchmark_realdata.csv", check.names = FALSE)
M <- tapply(b$uno_c, list(b$dataset, b$selector), identity)
sds$legacy_benchmark <- cbind(source = "legacy_benchmark", pair_sds(M))
cf <- list.files("results/realdata_cardinality", pattern = "^[A-Z]+\\.csv$", full.names = TRUE)
if (length(cf)) {
  cd <- do.call(rbind, lapply(cf, utils::read.csv, check.names = FALSE))
  for (K in sort(unique(cd$K))) {
    s <- cd[cd$K == K, ]
    MK <- tapply(s$uno_c, list(s$cohort, s$base), identity)
    if (ncol(MK) >= 2L) sds[[paste0("legacy_K", K)]] <- cbind(source = paste0("legacy_cardinality_K", K), pair_sds(MK))
  }
}
S <- do.call(rbind, sds)
RANKERS <- c("variance", "carss", "mim_surv")
pairs_of <- function(src, focus, rivals) S$source == src &
  ((S$a == focus & S$b %in% rivals) | (S$b == focus & S$a %in% rivals))
S$family <- ifelse(startsWith(S$source, "rcv_"), "other", "legacy")
S$family[pairs_of("rcv_main", "survfcbf_topk", RANKERS)] <- "headline"
for (K in c(5L, 10L, 20L))
  S$family[pairs_of(paste0("rcv_K", K), paste0("survfcbf_q2_K", K), paste0(RANKERS, "_K", K))] <- "matched"
k_sets <- sum(grepl("^K", names(RCV)))
stopifnot("headline needs 3 pairs" = !"main" %in% names(RCV) || sum(S$family == "headline") == 3L,
          "matched needs 3 pairs per K" = sum(S$family == "matched") == 3L * k_sets)
utils::write.csv(S, file.path(OUT, "pairwise_sd.csv"), row.names = FALSE)

Z80 <- stats::qnorm(0.975) + stats::qnorm(0.80); Z90 <- stats::qnorm(0.975) + stats::qnorm(0.90)
FAM <- intersect(c("headline", "matched"), unique(S$family))
sig <- do.call(rbind, lapply(FAM, function(f) {
  s <- S$sd_diff[S$family == f]
  data.frame(family = f, n_pairs = length(s), sigma = stats::median(s), sd_min = min(s), sd_max = max(s),
             delta_80_n11 = Z80 * stats::median(s) / sqrt(11), delta_90_n11 = Z90 * stats::median(s) / sqrt(11))
}))
utils::write.csv(sig, file.path(OUT, "sigma.csv"), row.names = FALSE)
cat("=== between-cohort SD of paired truncated Uno's C differences, 10-split means ===\n")
print(sig, digits = 3, row.names = FALSE)
print(S[S$family %in% FAM, c("family", "source", "a", "b", "sd_diff")], digits = 3, row.names = FALSE)
cat(sprintf("  legacy single-split untruncated median over all pairs, for context: %.3f\n",
            stats::median(S$sd_diff[S$family == "legacy"], na.rm = TRUE)))

## ---- 2. analytic n for a paired t-test ----------------------------------------
## optimistic and pessimistic use the smallest and largest SD of the family's pairs
n_needed <- function(delta, sigma, power = 0.8, alpha = 0.05)
  ceiling((stats::qnorm(1 - alpha / 2) + stats::qnorm(power))^2 * sigma^2 / delta^2)
deltas <- c(0.01, 0.02, 0.03, 0.05, 0.10)
an <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) data.frame(
  family = sig$family[i], sigma = sig$sigma[i], delta = deltas,
  n_80 = vapply(deltas, n_needed, numeric(1), sigma = sig$sigma[i], power = 0.80),
  n_90 = vapply(deltas, n_needed, numeric(1), sigma = sig$sigma[i], power = 0.90),
  n_80_optimistic = vapply(deltas, n_needed, numeric(1), sigma = sig$sd_min[i], power = 0.80),
  n_80_pessimistic = vapply(deltas, n_needed, numeric(1), sigma = sig$sd_max[i], power = 0.80))))
utils::write.csv(an, file.path(OUT, "analytic_n.csv"), row.names = FALSE)
cat("\n=== datasets needed, paired t-test (sigma = family median) ===\n")
print(an, digits = 3, row.names = FALSE)

## ---- 3. power of the Bayesian 0.95 rule under the argmax summary -----------------
## The posterior is summarised as in the Benavoli et al. reference implementation, by the
## fraction of draws in which each region holds the largest mass (summary = "argmax").
## Under that convention P_x converges to 1 for a real direction and the 0.95 rule behaves
## like an ordinary test, so this section reports its power: how often P_x > 0.95 fires,
## by N and delta, at the empirical sigma. The mean mass per region is not used here: under
## it P_x converges to P((Z1+Z2)/2 > rope), well below 0.95, however many datasets are
## added, and p_rival does not go to 0 for a real direction (0.099 at delta = 0.03).
## rate_rival_out, the share of benchmarks where P_rival < 0.05, behaves monotonically
## under argmax.
Ns    <- c(11L, 30L, 100L)
d_sim <- c(0.03, 0.05, 0.08, 0.12)
## 500 replicates per cell, so the Monte Carlo SE of a rate is at most 0.022
REPS  <- as.integer(Sys.getenv("POWER_REPS", "500")); NSIM <- 1500L
cat(sprintf("\n=== power of the argmax P>0.95 rule (ROPE +/-%.2f) ===\n", ROPE))
grid <- expand.grid(N = Ns, delta = d_sim, family = sig$family, stringsAsFactors = FALSE)
grid$sigma <- sig$sigma[match(grid$family, sig$family)]
grid$mean_px <- grid$rate_p95 <- grid$rate_rival_out <- NA_real_
for (i in seq_len(nrow(grid))) {
  px <- numeric(REPS); w95 <- 0L; wro <- 0L
  for (r in seq_len(REPS)) {
    z <- stats::rnorm(grid$N[i], mean = grid$delta[i], sd = grid$sigma[i])
    p <- bayes_signed_rank(z, rep(0, length(z)), rope = ROPE, nsim = NSIM, summary = "argmax")
    px[r] <- p[["p_x"]]
    if (!is.na(px[r]) && px[r] > 0.95) w95 <- w95 + 1L
    if (!is.na(p[["p_y"]]) && p[["p_y"]] < 0.05) wro <- wro + 1L
  }
  grid$mean_px[i] <- mean(px, na.rm = TRUE)
  grid$rate_p95[i] <- w95 / REPS
  grid$rate_rival_out[i] <- wro / REPS
}
grid$rate_p95_mcse <- sqrt(grid$rate_p95 * (1 - grid$rate_p95) / REPS)
grid$reps <- REPS
utils::write.csv(grid, file.path(OUT, "simulated_power.csv"), row.names = FALSE)
for (f in sig$family) for (dd in d_sim) {
  z <- grid[grid$family == f & grid$delta == dd, ]
  cat(sprintf("  %-8s sigma %.3f delta=%.2f  mean P_x: %s | rate(P_x>0.95): %s\n", f, z$sigma[1], dd,
              paste(sprintf("N=%d %.2f", z$N, z$mean_px), collapse = "  "),
              paste(sprintf("%.2f", z$rate_p95), collapse = "/")))
}
cat("\nReading: under the reference summary the 0.95 rule gains power with datasets, as a\n")
cat("decision rule should. The question for a benchmark is therefore the ordinary one --\n")
cat("whether N datasets give enough power at the gap of interest -- answered analytically\n")
cat("in section 2 above (Boulesteix et al. 2017 give the same calculation).\n")

saveRDS(list(sd = S, analytic = an, power = grid, sigma = sig), file.path(OUT, "power_analysis.rds"))
cat(sprintf("\ntables -> %s\nPOWER_DONE\n", OUT))
