# Canonical figures for the Application Note, derived from the result files. Every
# number quoted in the text, captions or abstract must come from here, so a corrected
# analysis propagates instead of leaving a stale value in the prose.
#
# Usage (from the repository root):  Rscript --no-init-file paper/make_statistics.R
# Output: paper/results/paper_statistics.json
suppressWarnings(suppressMessages(library(jsonlite)))
root <- Sys.getenv("SURVFCBF_ROOT", unset = ".")
setwd(root)
RES <- "experiments/benchmark/results"
OUT <- "paper/results"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## keep names: as.numeric() drops them, which silently turns a named selector-to-value map
## into an anonymous array in the JSON and makes the values unaddressable downstream
rd3 <- function(x) { v <- round(as.numeric(x), 3); names(v) <- names(x); v }

stats <- list(generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

## ---- cohorts -------------------------------------------------------------------
rd <- utils::read.csv(file.path(RES, "benchmark_realdata.csv"), check.names = FALSE)
coh <- unique(rd[, c("dataset", "n", "events", "p")])
stats$cohorts <- list(
  n_cohorts = nrow(coh), n_range = range(coh$n), events_range = range(coh$events),
  p_range = range(coh$p), total_patients = sum(coh$n), total_events = sum(coh$events))

## ---- main benchmark ------------------------------------------------------------
agg <- function(col) { v <- tapply(rd[[col]], rd$selector, function(x) mean(x, na.rm = TRUE))
                       as.list(rd3(v[order(-v)])) }
stats$benchmark <- list(
  uno_c = agg("uno_c"), stability = agg("stability"),
  parsimony = as.list(round(tapply(rd$parsimony, rd$selector, mean, na.rm = TRUE), 1)),
  failed_folds = as.list(round(tapply(rd$n_failed, rd$selector, mean, na.rm = TRUE), 2)),
  na_cohorts = as.list(tapply(rd$uno_c, rd$selector, function(x) sum(is.na(x)))),
  null_fit_folds = as.list(tapply(rd$n_null_fit, rd$selector, sum, na.rm = TRUE)))

## Bayesian signed-rank of survfcbf_topk against each rival, argmax summary, as written by
## realdata/analyze_cross_dataset.R (ROPE 0.02 on Uno's C, 0.05 on stability)
bsr <- utils::read.csv(file.path(RES, "realdata_crossdataset/bayes_signed_rank.csv"), check.names = FALSE)
stats$benchmark$bayes_vs_topk <- lapply(split(bsr, bsr$rival), function(z) list(
  uno_c = list(p_topk = z$C_topk, p_rope = z$C_rope, p_rival = z$C_rival, n = z$C_n,
               p_rival_mean_summary = z$C_rival_mean),
  stability = list(p_topk = z$S_topk, p_rope = z$S_rope, p_rival = z$S_rival)))

## 95% t interval of the paired Uno's C gap, top-K minus each ranker at its own cardinality,
## over the TCGA cohorts: the frequentist check of the headline that does not depend on how
## the Bayesian posterior is summarised
Mc <- tapply(rd$uno_c, list(rd$dataset, rd$selector), identity)
stats$benchmark$gap_ci_vs_topk <- lapply(stats::setNames(nm = c("variance", "carss", "mim_surv")), function(r) {
  d <- Mc[, "survfcbf_topk"] - Mc[, r]; d <- d[is.finite(d)]
  rd3(c(diff = mean(d), stats::t.test(d)$conf.int)) |> stats::setNames(c("diff", "lo", "hi")) |> as.list()
})

## ---- repeated cross-validation, the primary analysis -----------------------------------
## Uno's C truncated at the 80th percentile of each training fold's follow-up
## (realdata/uno_truncated.R), averaged over ten five-fold splits per cohort (seeds 1 to 10,
## seed 1 the published split). Headline: top-K against the rankers at the ten-percent
## cut-off (realdata/repeated_cv_main.R). Matched size: top-K at q = 2 against every rival
## cut to K = 5, 10, 20 genes (realdata/repeated_cv.R). Inference by
## realdata/repeated_cv_inference.R (argmax summary, ROPE 0.02 on C, 0.05 on stability).
rcv_read <- function(dir) do.call(rbind, lapply(list.files(file.path(RES, dir), pattern = "^[A-Z]+_rep[0-9]+\\.csv$",
                                                           full.names = TRUE), utils::read.csv, check.names = FALSE))
rcv_inf <- function(dir) {
  I <- utils::read.csv(file.path(RES, dir, "inference.csv"), check.names = FALSE)
  I$rival <- sub("_K[0-9]+$", "", I$rival)
  lapply(split(I, I$metric), function(z) lapply(split(z, z$rival), function(v) list(
    focus = rd3(v$focus_mean), rival = rd3(v$rival_mean), diff = rd3(v$mean_diff), se = rd3(v$se_diff),
    cohorts = paste0(v$wins_focus, "-", v$wins_rival), p_focus = rd3(v$p_focus), p_rope = rd3(v$p_rope),
    p_rival = rd3(v$p_rival), verdict = v$verdict, splits_same_verdict = v$reps_decided_same,
    sd_between_splits = rd3(v$sd_between_splits), sd_between_cohorts = rd3(v$sd_between_cohorts))))
}
rcv_means <- function(r) {
  r$sel <- sub("_K[0-9]+$", "", r$selector)
  cm <- aggregate(cbind(uno_c_trunc, uno_c, harrell_c, ibs, stability, parsimony) ~ cohort + sel, r, mean)
  lapply(stats::setNames(nm = c("uno_c_trunc", "uno_c", "ibs", "stability", "parsimony")), function(m)
    as.list(rd3(tapply(cm[[m]], cm$sel, mean))))
}
rcv_check <- function(r) stopifnot("repeated CV must have 11 cohorts x 10 splits and no failed fold" =
  nrow(unique(r[, c("cohort", "repetition")])) == 110L && all(r$n_failed == 0))
rm_ <- rcv_read("realdata_repeated_cv_main"); rcv_check(rm_)
cmn <- aggregate(uno_c_trunc ~ cohort + selector, rm_, mean)
Mt <- tapply(cmn$uno_c_trunc, list(cmn$cohort, cmn$selector), identity)
RANKERS <- c("variance", "carss", "mim_surv")
gap_ci <- lapply(stats::setNames(nm = RANKERS), function(r) { d <- Mt[, "survfcbf_topk"] - Mt[, r]
  as.list(stats::setNames(rd3(c(mean(d), stats::t.test(d)$conf.int)), c("diff", "lo", "hi"))) })
## integrated Brier score: no ROPE is defined for it, so it is reported descriptively, with
## its 95% t interval (positive = top-K worse)
cmi <- aggregate(ibs ~ cohort + selector, rm_, mean)
Mi <- tapply(cmi$ibs, list(cmi$cohort, cmi$selector), identity)
gap_ci_ibs <- lapply(stats::setNames(nm = RANKERS), function(r) { d <- Mi[, "survfcbf_topk"] - Mi[, r]
  list(diff = rd3(mean(d)), lo = rd3(stats::t.test(d)$conf.int[1]), hi = rd3(stats::t.test(d)$conf.int[2]),
       cohorts_topk_worse = sum(d > 0)) })
hi <- rcv_inf("realdata_repeated_cv_main")
stats$repeated_cv <- list(
  n_splits = length(unique(rm_$repetition)), truncation_quantile = 0.8,
  headline = list(means = rcv_means(rm_), inference = hi, gap_ci_trunc = gap_ci, gap_ci_ibs = gap_ci_ibs,
                  gap_range_trunc = rd3(range(vapply(hi$uno_c_trunc[RANKERS], function(v) v$diff, 0))),
                  null_fit_folds = as.list(tapply(rm_$n_null_fit, rm_$selector, sum))))
MK <- c(K5 = "realdata_repeated_cv_K5", K10 = "realdata_repeated_cv", K20 = "realdata_repeated_cv_K20")
mk <- lapply(MK, function(dir) { r <- rcv_read(dir); rcv_check(r)
  r$sel <- sub("_K[0-9]+$", "", r$selector)
  list(means = rcv_means(r), inference = rcv_inf(dir),
       null_fit_share = as.list(rd3(tapply(r$n_null_fit, r$sel, sum) / (5 * tapply(r$n_null_fit, r$sel, length))))) })
diff_range <- function(rival, metric = "uno_c_trunc") rd3(range(vapply(mk, function(z) z$inference[[metric]][[rival]]$diff, 0)))
verdicts <- function(rival, metric = "uno_c_trunc") unname(vapply(mk, function(z) z$inference[[metric]][[rival]]$verdict, ""))
RIVALS <- c(RANKERS, "mrmr", "mim_surv_q2", "mrmr_q2")
stats$repeated_cv$matched <- c(mk, list(
  range_trunc = lapply(stats::setNames(nm = RIVALS), diff_range),
  verdicts_trunc = lapply(stats::setNames(nm = RIVALS), verdicts),
  range_stability = lapply(stats::setNames(nm = RIVALS), diff_range, metric = "stability"),
  verdicts_stability = lapply(stats::setNames(nm = RIVALS), verdicts, metric = "stability")))

## sensitivity of the verdicts to the ROPE width and to the posterior summary
rope_sens <- function(dir) {
  z <- utils::read.csv(file.path(RES, dir, "rope_sensitivity.csv"), check.names = FALSE)
  z$rival <- sub("_K[0-9]+$", "", z$rival)
  lapply(split(z, z$summary), function(s) lapply(split(s, s$metric), function(m) lapply(split(m, m$rival), function(r)
    stats::setNames(as.list(r$verdict), paste0("rope_", r$rope)))))
}
stats$repeated_cv$headline$rope_sensitivity <- rope_sens("realdata_repeated_cv_main")
stats$repeated_cv$headline$p_rival_mean_summary <- local({
  z <- utils::read.csv(file.path(RES, "realdata_repeated_cv_main/rope_sensitivity.csv"), check.names = FALSE)
  z <- z[z$summary == "mean" & z$metric == "uno_c_trunc" & z$rope == 0.02, ]
  as.list(stats::setNames(rd3(z$p_rival), z$rival)) })
stats$repeated_cv$matched$rope_sensitivity <- lapply(MK, rope_sens)

## METABRIC, ten splits run locally (realdata/metabric_repeated_cv.R). Split 1 must
## reproduce the published benchmark's selections (stability); its C differs by the glmnet
## version (CESGA 5.0 for the published run, 4.1-10 here)
mrc <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_metabric/repeated_cv"), pattern = "^rep[0-9]+\\.csv$",
                                        full.names = TRUE), utils::read.csv, check.names = FALSE))
stopifnot("METABRIC repeated CV must have 10 splits, no failed fold, no null fit" =
            length(unique(mrc$repetition)) == 10L && all(mrc$n_failed == 0) && all(mrc$n_null_fit == 0))
mw <- function(m) tapply(mrc[[m]], list(mrc$repetition, mrc$selector), identity)
mgap <- function(m) lapply(stats::setNames(nm = RANKERS), function(r) { d <- mw(m)[, "survfcbf_topk"] - mw(m)[, r]
  list(diff = rd3(mean(d)), sd = rd3(stats::sd(d)), min = rd3(min(d)), max = rd3(max(d)),
       splits_topk_worse = if (m == "ibs") sum(d > 0) else sum(d < 0)) })
mb1 <- utils::read.csv(file.path(RES, "realdata_metabric/benchmark_metabric.csv"), check.names = FALSE)
m1 <- merge(mrc[mrc$repetition == 1L, ], mb1, by = "selector")
stopifnot("METABRIC split 1 no longer reproduces the published selections" =
            nrow(m1) == 4L && max(abs(m1$stability.x - m1$stability.y)) < 1e-9)
stats$repeated_cv$metabric <- list(
  means = lapply(stats::setNames(nm = c("uno_c_trunc", "uno_c", "ibs", "stability")), function(m)
    as.list(rd3(tapply(mrc[[m]], mrc$selector, mean)))),
  gap_trunc = mgap("uno_c_trunc"), gap_untrunc = mgap("uno_c"), gap_ibs = mgap("ibs"),
  split1_max_abs_c_diff_vs_published = signif(max(abs(m1$uno_c.x - m1$uno_c.y)), 2))

## untruncated Uno's C per fold (panel_redundancy.R folds): how often it is off by > 0.10
fo <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_redundancy"), pattern = "^[A-Z]+_folds\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
stats$repeated_cv$untruncated_off <- list(
  folds = nrow(fo), share_gap_over_0.10 = rd3(mean(abs(fo$uno_c - fo$uno_c_trunc) > 0.10)),
  share_gap_over_0.20 = rd3(mean(abs(fo$uno_c - fo$uno_c_trunc) > 0.20)))

## ---- power --------------------------------------------------------------------------
## sigma per family of comparisons (realdata/power_analysis.R): "headline" = top-K against
## the three rankers at their cut-off, "matched" = top-K at q = 2 against the same rankers
## at K = 5, 10, 20, both on the ten-split means of the truncated C. "legacy" rows keep the
## single-split, untruncated SDs of the earlier analysis for comparison only.
pw <- utils::read.csv(file.path(RES, "realdata_power/pairwise_sd.csv"), check.names = FALSE)
sg <- utils::read.csv(file.path(RES, "realdata_power/sigma.csv"), check.names = FALSE)
an <- utils::read.csv(file.path(RES, "realdata_power/analytic_n.csv"), check.names = FALSE)
sim <- utils::read.csv(file.path(RES, "realdata_power/simulated_power.csv"), check.names = FALSE)
stopifnot("power families missing" = setequal(sg$family, c("headline", "matched")))
by_fam <- function(z, col) setNames(as.list(z[[col]]), paste0("delta_", z$delta))
stats$power <- list(
  sigma = lapply(split(sg, sg$family), function(z) list(
    sigma = rd3(z$sigma), sd_min = rd3(z$sd_min), sd_max = rd3(z$sd_max), n_pairs = z$n_pairs,
    delta_80_n11 = rd3(z$delta_80_n11), delta_90_n11 = rd3(z$delta_90_n11))),
  datasets_needed_80 = lapply(split(an, an$family), by_fam, col = "n_80"),
  datasets_needed_90 = lapply(split(an, an$family), by_fam, col = "n_90"),
  datasets_needed_80_pessimistic = lapply(split(an, an$family), by_fam, col = "n_80_pessimistic"),
  ## power of the argmax P>0.95 rule at N=11 cohorts
  bayes_p95_power_n11 = lapply(split(sim[sim$N == 11, ], sim$family[sim$N == 11]), function(z)
    setNames(as.list(rd3(z$rate_p95)), paste0("delta_", z$delta))),
  reps = unique(sim$reps),
  legacy_sd_median = rd3(stats::median(pw$sd_diff[pw$family == "legacy"])))

## ---- matched cardinality -------------------------------------------------------
cf <- list.files(file.path(RES, "realdata_cardinality"), pattern = "^[A-Z]+\\.csv$", full.names = TRUE)
cd <- do.call(rbind, lapply(cf, utils::read.csv, check.names = FALSE))
ca <- aggregate(uno_c ~ base + K, cd, function(x) mean(x, na.rm = TRUE))
cs <- aggregate(stability ~ base + K, cd, function(x) mean(x, na.rm = TRUE))
by_K <- function(a, col) lapply(split(a, a$base), function(z) setNames(as.list(rd3(z[[col]])), paste0("K", z$K)))
inf <- utils::read.csv(file.path(RES, "realdata_cardinality/cardinality_inference.csv"), check.names = FALSE)
vs <- function(rival) { z <- inf[inf$rival == rival, ]
  lapply(split(z, z$K), function(w) list(diff = rd3(w$mean_diff), se = rd3(w$se_diff),
                                         cohorts = paste0(w$wins_focus, "-", w$wins_rival),
                                         p_rival = rd3(w$p_rival), verdict = w$verdict)) }
stats$cardinality <- list(
  uno_c_by_K = by_K(ca, "uno_c"), stability_by_K = by_K(cs, "stability"),
  vs_variance = vs("variance"), vs_mim_surv = vs("mim_surv"), vs_carss = vs("carss"))
## top-K returns fewer than K genes once the cohort runs out of predominant features, so
## only the sizes where survfcbf's mean panel reaches 90% of K count as matched; the note
## quotes the range of the paired gaps over those sizes
fk <- aggregate(parsimony ~ K, cd[cd$base == "survfcbf", ], mean)
mK <- fk$K[fk$parsimony >= 0.9 * fk$K]
rng <- function(rival) rd3(range(inf$mean_diff[inf$rival == rival & inf$K %in% mK]))
stats$cardinality$matched <- list(
  K = mK, survfcbf_panel_by_K = setNames(as.list(round(fk$parsimony, 1)), paste0("K", fk$K)),
  range_vs_variance = rng("variance"), range_vs_mim_surv = rng("mim_surv"),
  range_vs_carss = rng("carss"),
  all_undecided = all(inf$verdict[inf$K %in% mK] == "undecided at n=11"))

## ---- matched size: redundancy, mRMR, censoring and q (realdata/panel_redundancy.R) ------
RR <- file.path(RES, "realdata_redundancy")
rr <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+\\.csv$", full.names = TRUE),
                            utils::read.csv, check.names = FALSE))
rv <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+_reverse\\.csv$", full.names = TRUE),
                            utils::read.csv, check.names = FALSE))
rr$K <- as.integer(sub(".*_K", "", rr$selector)); rr$family <- sub("_K[0-9]+$", "", rr$selector)
## same folds, seed and learner as the cardinality curve, so the shared selectors must
## reproduce it exactly
cdf <- cd; cdf$family <- ifelse(cdf$base == "survfcbf", "survfcbf_q2", cdf$base)
cmp <- merge(rr, cdf, by = c("cohort", "family", "K"))
stopifnot("panel_redundancy.R no longer reproduces the cardinality curve" =
            nrow(cmp) == 4L * 3L * length(unique(rr$cohort)) &&
            max(abs(cmp$uno_c.x - cmp$uno_c.y)) < 1e-9 &&
            max(abs(cmp$stability.x - cmp$stability.y)) < 1e-9)
k10 <- rr[rr$K == 10L, ]
fam <- function(z, col) as.list(rd3(tapply(z[[col]], z$family, mean, na.rm = TRUE)))
at10 <- function(f, col) { z <- k10[k10$family == f, ]; z[[col]][order(z$cohort)] }
mi <- utils::read.csv(file.path(RR, "mrmr_inference.csv"), check.names = FALSE)
stats$redundancy <- list(
  n_cohorts = length(unique(rr$cohort)),
  mean_abs_r_K10 = fam(k10, "mean_abs_r"), share_r05_K10 = fam(k10, "share_r05"),
  uno_c_K10 = fam(k10, "uno_c"), stability_K10 = fam(k10, "stability"),
  cohorts_survfcbf_below_mim = sum(at10("survfcbf_q2", "mean_abs_r") < at10("mim_surv", "mean_abs_r")),
  cohorts_survfcbf_below_mim_q2 = sum(at10("survfcbf_q2", "mean_abs_r") < at10("mim_surv_q2", "mean_abs_r")),
  cohorts_survfcbf_below_mrmr_q2 = sum(at10("survfcbf_q2", "mean_abs_r") < at10("mrmr_q2", "mean_abs_r")),
  failed_folds = sum(rr$n_failed),
  vs_mrmr = lapply(split(mi, mi$metric), function(z) lapply(split(z, z$K), function(v) list(
    diff = rd3(v$mean_diff), se = rd3(v$se_diff), cohorts = paste0(v$wins_focus, "-", v$wins_rival),
    p_focus = rd3(v$p_focus), p_rival = rd3(v$p_rival), verdict = v$verdict))))
rvm <- tapply(rv$harrell_c_censoring, rv$selector, mean)
## each panel's Cox model on censoring against 100 random ten-gene panels of its cohort
rnd <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+_reverse_random\\.csv$", full.names = TRUE),
                             utils::read.csv, check.names = FALSE))
stopifnot("expected 100 random panels per cohort" = all(table(rnd$cohort) == 100L))
rv$pct <- mapply(function(co, v) mean(rnd$harrell_c_censoring[rnd$cohort == co] < v), rv$cohort, rv$harrell_c_censoring)
rvs <- rv[rv$selector != "random_K10", ]
q2p <- rvs[rvs$selector == "survfcbf_q2_K10", ]
mcc <- utils::read.csv(file.path(RES, "realdata_metabric/censoring_cox.csv"), check.names = FALSE)
stats$censoring_cox <- list(harrell_c = as.list(rd3(rvm)), failed_folds = sum(rv$n_failed),
                            null_fit_folds = sum(rv$n_null_fit),
                            random_panels = 100L, random_mean = rd3(mean(rnd$harrell_c_censoring)),
                            cohorts_above_random_p95 = as.list(tapply(rvs$pct > 0.95, rvs$selector, sum)),
                            survfcbf_q2_max_percentile = rd3(max(q2p$pct)),
                            survfcbf_q2_max_percentile_cohort = q2p$cohort[which.max(q2p$pct)],
                            metabric_percentile = as.list(stats::setNames(rd3(mcc$percentile_among_random), mcc$selector)),
                            metabric_harrell_c = as.list(stats::setNames(rd3(mcc$harrell_c_censoring), mcc$selector)))

## ---- stability threshold (realdata/stability_theta.R): theta 0.9 must reproduce the panels
th <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_stability_theta"), pattern = "\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
stopifnot("theta = 0.9 no longer reproduces panel_redundancy.R" =
            nrow(th) == 3L * 77L && max(abs(th$stability[th$theta == 0.9] - th$stability_main_run[th$theta == 0.9])) < 1e-12)
th$family <- sub("_K10$", "", th$selector)
thw <- function(rival) lapply(split(th, th$theta), function(z) {
  a <- z[z$family == "survfcbf_q2", ]; b <- z[z$family == rival, ]; m <- merge(a, b, by = "cohort")
  list(diff = rd3(mean(m$stability.x - m$stability.y)), cohorts_lower = sum(m$stability.x < m$stability.y)) })
stats$stability_theta <- list(
  thetas = sort(unique(th$theta)), draws = as.list(stats::setNames(tapply(th$draws, th$theta, unique), sort(unique(th$theta)))),
  mean = lapply(split(th, th$theta), function(z) as.list(rd3(tapply(z$stability, z$family, mean)))),
  vs_mim_surv_q2 = thw("mim_surv_q2"), vs_mrmr_q2 = thw("mrmr_q2"), vs_mrmr = thw("mrmr"))

## ---- survival target: martingale against IPCW-binned (realdata/target_sensitivity.R) -----
ts <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_target_sensitivity"), pattern = "^[A-Z]+\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
stopifnot("target sensitivity: 11 cohorts x 2 targets, no failed fold" = nrow(ts) == 22L && all(ts$n_failed == 0))
tw <- function(col) tapply(ts[[col]], list(ts$cohort, ts$target), identity)
stats$target <- list(
  means = lapply(stats::setNames(nm = c("uno_c_trunc", "uno_c", "stability", "mean_abs_r", "harrell_c_censoring")),
                 function(m) as.list(rd3(tapply(ts[[m]], ts$target, mean)))),
  diff_trunc_ipcw_minus_martingale = rd3(mean(tw("uno_c_trunc")[, "ipcw"] - tw("uno_c_trunc")[, "martingale"])),
  jaccard_mean = rd3(mean(ts$jaccard_between_targets)), jaccard_max = rd3(max(ts$jaccard_between_targets)),
  ## the martingale gradient: does the panel's censoring C rise with the censored share of
  ## the upper class? Spearman over cohorts, per target
  gradient_spearman = lapply(split(ts, ts$target), function(z)
    rd3(stats::cor(z$share_censored_upper, z$harrell_c_censoring, method = "spearman"))),
  brca_censoring_c = as.list(rd3(stats::setNames(ts$harrell_c_censoring[ts$cohort == "BRCA"], ts$target[ts$cohort == "BRCA"]))))

## ---- centre and follow-up (realdata/center_followup.R) -------------------------------
ce <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_center"), pattern = "^[A-Z]+\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
stats$centre <- list(
  n_centres_range = range(ce$n_centres), largest_share_range = rd3(range(ce$largest_centre_share)),
  c_centre_censoring = rd3(mean(ce$c_centre_censoring)), c_centre_censoring_range = rd3(range(ce$c_centre_censoring)),
  c_centre_death = rd3(mean(ce$c_centre_death)), c_centre_death_range = rd3(range(ce$c_centre_death)),
  cohorts_followup_differs_p001 = sum(ce$kruskal_p_followup_censored < 0.001),
  median_followup_days_range = range(ce$median_followup_reverse_km),
  time_max_days_range = range(ce$time_max))

## ---- METABRIC clinical reference (realdata/metabric_clinical.R) ----------------------
cl <- utils::read.csv(file.path(RES, "realdata_metabric/clinical_reference.csv"), check.names = FALSE)
clr <- function(m) cl[cl$model == m, ]
stats$metabric_clinical <- list(
  n = clr("clinical")$n, events = clr("clinical")$events,
  uno_c_trunc = as.list(stats::setNames(rd3(cl$uno_c_trunc[cl$model != "clinical_plus_panel_minus_clinical"]),
                                        cl$model[cl$model != "clinical_plus_panel_minus_clinical"])),
  added_by_panel_trunc = rd3(clr("clinical_plus_panel_minus_clinical")$uno_c_trunc),
  added_by_panel_trunc_se = rd3(clr("clinical_plus_panel_minus_clinical")$uno_c_trunc_se),
  added_by_panel_untrunc = rd3(clr("clinical_plus_panel_minus_clinical")$uno_c),
  added_by_panel_ibs = rd3(clr("clinical_plus_panel_minus_clinical")$ibs))

## ---- simulation at matched size (sim_matched_size.R, the 300 datasets of the grid) -------
sm <- utils::read.csv(file.path(RES, "raw_sim_matched_size.csv"), check.names = FALSE)
sm$K <- as.integer(sub(".*_K", "", sm$selector)); sm$family <- sub("_K[0-9]+$", "", sm$selector)
stopifnot("matched-size simulation lost folds" = all(sm$n_failed == 0))
sim_fam <- function(col) lapply(split(sm, sm$K), function(z) as.list(rd3(tapply(z[[col]], z$family, mean))))
## the note quotes recovery to two decimals: round the raw means once, never a rounded value
## again (0.3447 rounded to 0.345 would print as 0.35)
## best rival at its default binning (the rankers and mRMR at q = 10); MIM and mRMR at
## q = 2, top-K's own binning, are reported apart because mRMR at q = 2 is top-K's twin
best_rival <- function(col, d = 3) lapply(split(sm, sm$K), function(z) {
  v <- tapply(z[[col]], z$family, mean); round(max(v[names(v) != "survfcbf_q2" & !grepl("_q2$", names(v))]), d) })
own <- function(col, d) lapply(split(sm, sm$K), function(z) round(mean(z[[col]][z$family == "survfcbf_q2"]), d))
stats$sim_matched <- list(
  n_datasets = nrow(unique(sm[, c("scenario", "replicate")])), K = sort(unique(sm$K)),
  sensitivity = sim_fam("sensitivity"), redundant_copies = sim_fam("redundant_copies"),
  fdr = sim_fam("fdr"), exact = sim_fam("exact"), uno_c = sim_fam("uno_c"),
  stability = sim_fam("stability"), best_rival_sensitivity = best_rival("sensitivity"),
  sensitivity_2dp = list(survfcbf = own("sensitivity", 2), best_rival = best_rival("sensitivity", 2),
                         mrmr_q2 = lapply(split(sm, sm$K), function(z) round(mean(z$sensitivity[z$family == "mrmr_q2"]), 2)),
                         mim_surv_q2 = lapply(split(sm, sm$K), function(z) round(mean(z$sensitivity[z$family == "mim_surv_q2"]), 2))),
  uno_c_trunc = sim_fam("uno_c_trunc"),
  survfcbf_redundant_copies_total = sum(sm$redundant_copies[sm$family == "survfcbf_q2"]),
  null_fit_folds = as.list(tapply(sm$n_null_fit, sm$family, sum)))

## ---- the evaluation traps ------------------------------------------------------
mv <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_modules_v3"),
                                       pattern = "^[A-Z]+\\.rds$", full.names = TRUE),
                            function(f) readRDS(f)$summary))
ml <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_modules_leaksafe"),
                                       pattern = "rds$", full.names = TRUE), readRDS))
mg <- merge(mv[, c("cohort", "stab_module_su", "stab_module_car")],
            ml[, c("cohort", "stab_module_su_ls", "stab_module_car_ls", "stab_gene")], by = "cohort")
stats$leakage <- list(
  module_su_fixed = rd3(mean(mg$stab_module_su)), module_su_leaksafe = rd3(mean(mg$stab_module_su_ls)),
  module_car_fixed = rd3(mean(mg$stab_module_car)), module_car_leaksafe = rd3(mean(mg$stab_module_car_ls)),
  gene_level = rd3(mean(mg$stab_gene)),
  cohorts_below_gene_level = sum(mg$stab_module_car_ls < mg$stab_gene), n_cohorts = nrow(mg),
  inflation_factor = rd3(mean(mg$stab_module_car) / mean(mg$stab_module_car_ls)))

ec <- utils::read.csv(file.path(RES, "realdata_metabric/event_curve.csv"), check.names = FALSE)
## benchmark_metabric.csv is the lambda.min nested benchmark on the complete profile
## (sha256 46bee94b..., CESGA run of 2026-10-30): the source for the cohort's n, p and
## the ranker panel size the note quotes, which until now were quoted without a file
mb <- utils::read.csv(file.path(RES, "realdata_metabric/benchmark_metabric.csv"), check.names = FALSE)
stopifnot("METABRIC benchmark has null-fit folds; investigate before quoting it" =
            all(mb$n_null_fit == 0))
stats$metabric <- list(
  n = unique(mb$n), events = unique(mb$events), p = unique(mb$p),
  ranker_panel = unique(mb$parsimony[mb$selector %in% c("variance", "carss", "mim_surv")]),
  topk_panel = mb$parsimony[mb$selector == "survfcbf_topk"],
  events_range = range(ec$events),
  survfcbf_stability_range = rd3(range(ec$`stability.survfcbf_topk`, na.rm = TRUE)),
  variance_stability_range = rd3(range(ec$stability.variance, na.rm = TRUE)))

## ---- mechanism and biology -----------------------------------------------------
ca2 <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_car_ablation"),
                                        pattern = "rds$", full.names = TRUE), readRDS))
stats$car_ablation <- list(
  su_mean = rd3(mean(ca2$stab_su)), car_mean = rd3(mean(ca2$stab_car)),
  cohorts_car_better = sum(ca2$stab_car > ca2$stab_su), n_cohorts = nrow(ca2),
  consensus_su = round(mean(ca2$consensus_su), 1), consensus_car = round(mean(ca2$consensus_car), 1))

gs <- readRDS(file.path(RES, "realdata_gsea/gsea.rds"))
gsea_hits <- vapply(sort(unique(gs$selector)), function(sl) {
  r <- gs[gs$selector == sl & !is.na(gs$padj), ]
  sum(tapply(r$padj < 0.05, r$cohort, any, na.rm = TRUE)) }, numeric(1))

## ORA at K=10, COMPUTED from the raw per-set tests (never hard-coded: an earlier
## version wrote these by hand and the number-checking script could not catch it).
## Both correction policies are reported, labelled, because they answer different
## questions: global BH treats every set x signature x cohort test as one family,
## per-signature BH (the ORA convention) corrects within each cohort x selector.
## The Application Note must name the policy it quotes.
ora_hits <- function(dir) {
  o <- utils::read.csv(file.path(RES, dir, "ora_full.csv"), check.names = FALSE)
  o <- o[o$K == 10L, ]
  o$q_sig <- stats::ave(o$p, paste(o$cohort, o$selector), FUN = function(p) stats::p.adjust(p, "BH"))
  sel <- sort(unique(o$selector))
  list(
    global_bh = as.list(vapply(sel, function(sl) { z <- o[o$selector == sl, ]
      sum(tapply(z$q < 0.05, z$cohort, any, na.rm = TRUE)) }, numeric(1))),
    per_signature_bh = as.list(vapply(sel, function(sl) { z <- o[o$selector == sl, ]
      sum(tapply(z$q_sig < 0.05, z$cohort, any, na.rm = TRUE)) }, numeric(1))))
}
stats$enrichment <- list(
  gsea_cohorts_enriched = as.list(gsea_hits),
  ora_topk_full_space = ora_hits("realdata_enrichment"),
  ora_topk_coding_only = ora_hits("realdata_enrichment_coding"))

rt <- utils::read.csv(file.path(RES, "realdata_runtime/runtime_scaling.csv"), check.names = FALSE)
fit <- stats::lm(log(survfcbf_topk) ~ log(p), data = rt)
stats$runtime <- list(scaling_exponent = rd3(unname(coef(fit)[2])),
                      seconds_at_p20k = rd3(rt$survfcbf_topk[which.max(rt$p)]),
                      variance_seconds_at_p20k = rd3(rt$variance[which.max(rt$p)]))

## ---- simulation grid (cesga_core: 10 configurations x 30 replicates, outer = 5) -----
## Two different failures, kept apart because they were once conflated. A replicate has
## an EMPTY panel when no outer fold returned a feature (recovery is then NA, since
## benchmark_recovery() scores only non-empty folds, and parsimony is NA too). A fold
## FAILS when the panel has fewer than the two features the ridge Cox needs, so a
## one-gene panel fails without being empty. The App Note quotes the empty-panel share.
sg <- utils::read.csv(file.path(RES, "raw_cesga_core.csv"), check.names = FALSE)
OUTER <- 5L                                               # cesga_core tier in grid_lib.R
## The three one-at-a-time axes share the reference configuration (n = 600, 50% censoring,
## p = 2000). The events axis and the censoring axis both include it, as A_n600 and
## A_cens50, and the replicate seed depends only on the replicate, so the two labels hold
## the same 30 datasets. Count them once: the grid has 10 distinct configurations and 300
## distinct datasets, not 11 and 330.
DUP <- c(dup = "A_cens50", of = "A_n600")
same_rows <- function(a, b) {
  key <- setdiff(names(sg), c("scenario", "map", "runtime", "runtime_se"))
  a <- a[order(a$replicate, a$selector), key]; b <- b[order(b$replicate, b$selector), key]
  isTRUE(all.equal(a, b, check.attributes = FALSE))
}
stopifnot("A_cens50 is no longer a copy of A_n600: revisit the deduplication" =
            same_rows(sg[sg$scenario == DUP[["dup"]], ], sg[sg$scenario == DUP[["of"]], ]))
sg <- sg[sg$scenario != DUP[["dup"]], ]
stopifnot("empty panel: recovery and parsimony must agree" =
            identical(is.na(sg$sensitivity), is.na(sg$parsimony)))
sim_by <- function(f) as.list(vapply(c("survfcbf", "survfcbf_cvq"), function(sl)
  rd3(f(sg[sg$selector == sl, ])), numeric(1)))
stats$simulation <- list(
  n_scenarios = length(unique(sg$scenario)), n_replicates = sum(sg$selector == "survfcbf"),
  empty_panel_share = sim_by(function(z) mean(is.na(z$sensitivity))),
  all_folds_failed_share = sim_by(function(z) mean(z$n_failed >= OUTER)),
  failed_fold_share = sim_by(function(z) mean(z$n_failed / OUTER)))

## ---- correctness against the reference FCBF package ---------------------------------
xv <- utils::read.csv(file.path(RES, "validation/fcbf_crossvalidation.csv"), check.names = FALSE)
xv_row <- function(check) xv[xv$check == check, ]
topk_rows <- xv[grepl("^topk_", xv$check), ]
stats$validation <- list(
  su_pairs = xv_row("su_random_pairs")$comparisons,
  su_max_abs_diff = signif(xv_row("su_random_pairs")$max_abs_diff, 3),
  selections_identical = xv_row("selections_nondegenerate")$comparisons - xv_row("selections_nondegenerate")$mismatches,
  selections_mismatches = xv_row("selections_nondegenerate")$mismatches,
  selections_skipped_exact_ties = xv_row("selections_nondegenerate")$skipped_ties,
  tie_sweep_datasets = xv_row("selections_with_ties")$comparisons,
  tie_sweep_mismatches = xv_row("selections_with_ties")$mismatches,
  tie_sweep_mismatches_at_exact_tie = xv_row("selections_with_ties")$mismatches_at_tie,
  scdengue_genes = as.integer(sub(" genes.*", "", xv_row("selections_scdengue")$detail)),
  scdengue_thresholds = xv_row("selections_scdengue")$comparisons,
  scdengue_mismatches = xv_row("selections_scdengue")$mismatches,
  topk_comparisons = sum(topk_rows$comparisons), topk_mismatches = sum(topk_rows$mismatches))

jsonlite::write_json(stats, file.path(OUT, "paper_statistics.json"), pretty = TRUE, auto_unbox = TRUE)
cat("wrote", file.path(OUT, "paper_statistics.json"), "\n")
cat(sprintf("  cohorts %d | headline sigma %.3f | datasets for delta=0.03: %d | leakage inflation %.1fx\n",
            stats$cohorts$n_cohorts, stats$power$sigma$headline$sigma,
            stats$power$datasets_needed_80$headline$delta_0.03, stats$leakage$inflation_factor))
cat("STATISTICS_DONE\n")
