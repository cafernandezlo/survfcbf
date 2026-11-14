# Check that the figures quoted in the Application Note still match the computed results.
# Rewriting a sentence that contains a number, or rerunning an analysis, breaks the
# correspondence silently. This fails loudly instead.
#
# Run after ANY edit to the text or to the pipeline:
#   Rscript --no-init-file paper/make_statistics.R
#   Rscript --no-init-file paper/check_manuscript_numbers.R
#
# Each check names a claim, the value the results produce, and how it should appear in
# the text. A "rule" is a claim the text makes in words (e.g. "no comparison is decided")
# checked against the results directly. Exit status is non-zero if the note misses a
# claim, so this can gate a commit.
root <- Sys.getenv("SURVFCBF_ROOT", unset = "."); setwd(root)
suppressMessages(library(jsonlite))
s <- jsonlite::fromJSON("paper/results/paper_statistics.json")

## LaTeX writes 2{,}398 and 48\%, the results give 2398 and 0.482. A document can span
## several files (the note's text is split into pieces shared with the preprint).
read_tex <- function(path) {
  tex <- paste(unlist(lapply(path, readLines, warn = FALSE)), collapse = " ")
  gsub("\\\\%", "%", gsub("\\{,\\}", ",", tex))
}
fmt  <- function(x, d = 3) formatC(as.numeric(x), format = "f", digits = d)
big  <- function(x) formatC(round(as.numeric(x)), format = "d", big.mark = ",")
pct  <- function(x) paste0(round(100 * as.numeric(x)), "%")
WORD <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
          "eleven", "twelve")
word <- function(n) WORD[as.integer(n)]
num  <- function(name, value) list(name = name, value = value, type = "num")
## a count in words is only meaningful in its phrase ("ten" alone also matches "ten-gene")
wrd  <- function(name, n, phrase) list(name = name, value = sub("%s", word(n), phrase, fixed = TRUE),
                                       type = "word")
rule <- function(name, ok, shown) list(name = name, value = shown, type = "rule", ok = isTRUE(ok))

found <- function(ck, tex) switch(ck$type,
  ## match the value as a standalone token so 0.055 does not match 0.0551
  num  = grepl(paste0("(^|[^0-9.,])", gsub("\\.", "\\\\.", ck$value), "([^0-9]|$)"), tex),
  word = grepl(paste0("\\b", ck$value, "\\b"), tex, ignore.case = TRUE),
  rule = ck$ok)

## ---- Application Note ----------------------------------------------------------------
## The primary analysis is the repeated cross-validation (ten five-fold splits, Uno's C
## truncated at the 80th percentile of the training follow-up), s$repeated_cv. The panel
## sizes, failures, runtime, simulation of the permutation mode, enrichment and validation
## come from the published single split and do not depend on the metric.
rankers <- c("variance", "carss", "mim_surv")          # the three at their own cardinality
rc  <- s$repeated_cv
hl  <- rc$headline
hm  <- hl$means$uno_c_trunc
## the gap is the mean paired difference over cohorts, the number the supplement gives with
## its interval, not the difference of the two rounded means
gaps <- -vapply(hl$gap_ci_trunc[rankers], function(g) g$diff, numeric(1))
hv   <- function(r, m = "uno_c_trunc") hl$inference[[m]][[r]]$verdict
mtc  <- rc$matched
Ks   <- as.integer(sub("^K", "", grep("^K[0-9]+$", names(mtc), value = TRUE)))
mgap <- -vapply(rc$metabric$gap_trunc[rankers], function(g) g$diff, numeric(1))
pw   <- s$power
red  <- s$redundancy
k10s <- mtc$K10$means$stability
th   <- s$stability_theta$mean
cz   <- s$censoring_cox
sm2  <- s$sim_matched$sensitivity_2dp
appnote_checks <- list(
  wrd("TCGA cohorts",                        s$cohorts$n_cohorts, "%s TCGA"),
  num("TCGA features low",                   big(s$cohorts$p_range[1])),
  num("TCGA features high",                  big(s$cohorts$p_range[2])),
  num("TCGA patients low",                   as.character(s$cohorts$n_range[1])),
  num("TCGA patients high",                  as.character(s$cohorts$n_range[2])),
  num("TCGA events low",                     as.character(s$cohorts$events_range[1])),
  num("TCGA events high",                    as.character(s$cohorts$events_range[2])),
  num("METABRIC events",                     big(s$metabric$events_range[2])),
  num("panel size, rankers",                 big(s$benchmark$parsimony$variance)),
  num("panel size, permutation mode",        fmt(s$benchmark$parsimony$survfcbf, 1)),
  num("panel size, penalised Cox",           fmt(s$benchmark$parsimony$cox_glmnet_l1, 1)),
  wrd("cohorts lost, permutation mode",      s$benchmark$na_cohorts$survfcbf, "%s of eleven"),
  rule("top-K never fails",                  s$benchmark$failed_folds$survfcbf_topk == 0 &&
                                               s$benchmark$na_cohorts$survfcbf_topk == 0, "0 failed folds"),
  num("empty panels, simulation",            pct(s$simulation$empty_panel_share$survfcbf)),
  ## anchored to its phrase: a bare "300" also matches the validation "sweep of 300 datasets"
  num("simulated datasets",                  sprintf("%d simulated datasets", s$simulation$n_replicates)),
  num("runtime exponent",                    fmt(s$runtime$scaling_exponent, 2)),
  num("seconds at p = 20k",                  fmt(s$runtime$seconds_at_p20k, 1)),
  ## metric and design
  num("truncation percentile",               sprintf("%dth percentile", round(100 * rc$truncation_quantile))),
  wrd("repetitions of cross-validation",     rc$n_splits, "%s repetitions"),
  num("folds where untruncated C is off",    sprintf("%.1f%%", 100 * rc$untruncated_off$share_gap_over_0.10)),
  ## power, headline family
  num("SD of paired C differences",          fmt(pw$sigma$headline$sigma)),
  num("detectable difference, paired test",  fmt(pw$sigma$headline$delta_80_n11, 2)),
  num("power of the 0.95 rule at 0.05",      pct(pw$bayes_p95_power_n11$headline$delta_0.05)),
  rule("the 0.95 rule is weak at 0.03, near 80% only at 0.05",
       pw$bayes_p95_power_n11$headline$delta_0.03 < 0.5 && pw$bayes_p95_power_n11$headline$delta_0.05 > 0.75,
       sprintf("%.2f / %.2f", pw$bayes_p95_power_n11$headline$delta_0.03, pw$bayes_p95_power_n11$headline$delta_0.05)),
  ## headline comparison
  num("truncated C, top-K",                  fmt(hm$survfcbf_topk)),
  num("truncated C, lowest ranker",          fmt(min(unlist(hm[rankers])))),
  num("truncated C, highest ranker",         fmt(max(unlist(hm[rankers])))),
  num("gap to rankers, smallest",            fmt(min(gaps))),
  num("gap to rankers, largest",             fmt(max(gaps))),
  num("gap to rankers, about (abstract)",    fmt(mean(gaps), 2)),
  rule("95% intervals of the gap exclude zero", all(vapply(hl$gap_ci_trunc, function(g) g$hi < 0, logical(1))),
       "upper limits below zero"),
  rule("0.95 rule decides all three at ROPE 0.02", all(vapply(rankers, hv, "") == "rival better"), "3 of 3"),
  rule("none decided at ROPE 0.03",          all(vapply(rankers, function(r)
         hl$rope_sensitivity$argmax$uno_c_trunc[[r]]$rope_0.03, "") != "rival better"), "0 of 3"),
  ## METABRIC
  num("METABRIC gap, smallest",              fmt(min(mgap))),
  num("METABRIC gap, largest",               fmt(max(mgap))),
  rule("METABRIC gap in all ten splits",     all(vapply(rc$metabric$gap_trunc, function(g) g$splits_topk_worse == 10L, logical(1))),
       "10 of 10"),
  num("METABRIC clinical, panel adds",       fmt(s$metabric_clinical$added_by_panel_trunc)),
  num("METABRIC clinical, its SE",           fmt(s$metabric_clinical$added_by_panel_trunc_se)),
  ## matched panel size, ten splits
  num("matched sizes, largest K",            sprintf("between %d and %d genes", min(Ks), max(Ks))),
  num("matched gap to variance, smallest",   fmt(mtc$range_trunc$variance[1])),
  num("matched gap to variance, largest",    fmt(mtc$range_trunc$variance[2])),
  num("matched gap to MIM, smallest",        fmt(-mtc$range_trunc$mim_surv[2])),
  num("matched gap to MIM, largest",         fmt(-mtc$range_trunc$mim_surv[1])),
  rule("ahead of variance, behind MIM",      all(mtc$range_trunc$variance > 0) && all(mtc$range_trunc$mim_surv < 0),
       "signs of the matched gaps"),
  rule("neither variance nor MIM decided",   all(c(mtc$verdicts_trunc$variance, mtc$verdicts_trunc$mim_surv) == "undecided at n=11"),
       "6 undecided"),
  rule("practically equivalent to mRMR at every size", all(mtc$verdicts_trunc$mrmr == "practically equivalent"),
       paste(Ks, collapse = ",")),
  ## redundancy, on the same two-bin target
  num("in-panel |r| at K=10, survfcbf",      fmt(red$mean_abs_r_K10$survfcbf_q2)),
  num("in-panel |r| at K=10, MIM q=2",       fmt(red$mean_abs_r_K10$mim_surv_q2)),
  num("in-panel |r| at K=10, mRMR q=2",      fmt(red$mean_abs_r_K10$mrmr_q2)),
  num("in-panel |r| at K=10, CAR",           fmt(red$mean_abs_r_K10$carss)),
  rule("pruning halves the MIM redundancy",  abs(red$mean_abs_r_K10$survfcbf_q2 / red$mean_abs_r_K10$mim_surv_q2 - 0.5) < 0.1,
       sprintf("ratio %.2f", red$mean_abs_r_K10$survfcbf_q2 / red$mean_abs_r_K10$mim_surv_q2)),
  wrd("cohorts below MIM q=2",               red$cohorts_survfcbf_below_mim_q2, "all %s cohorts"),
  rule("mRMR and CAR at a similar or lower level", red$mean_abs_r_K10$mrmr_q2 < red$mean_abs_r_K10$survfcbf_q2 &&
         red$mean_abs_r_K10$carss < red$mean_abs_r_K10$survfcbf_q2, "both below survfcbf"),
  ## stability at ten genes, ten splits, and its threshold
  num("stability K=10, survfcbf",            fmt(k10s$survfcbf_q2)),
  num("stability K=10, MIM q=2",             fmt(k10s$mim_surv_q2)),
  num("stability K=10, mRMR q=2",            fmt(k10s$mrmr_q2)),
  num("stability K=10, mRMR",                fmt(k10s$mrmr)),
  num("stability K=10, CAR",                 fmt(k10s$carss)),
  rule("matches mRMR on the same target",    mtc$K10$inference$stability$mrmr_q2$verdict == "practically equivalent",
       mtc$K10$inference$stability$mrmr_q2$verdict),
  rule("ordering holds at every threshold",  all(vapply(th, function(z) z$survfcbf_q2 < z$mim_surv_q2 &&
         z$survfcbf_q2 > z$mrmr && z$survfcbf_q2 < z$carss, logical(1))) &&
         k10s$survfcbf_q2 < k10s$mim_surv_q2 && k10s$survfcbf_q2 > k10s$mrmr && k10s$survfcbf_q2 < k10s$carss,
       paste(names(th), collapse = ",")),
  num("lowest similarity threshold",         sprintf("from %s to %s", min(s$stability_theta$thetas), max(s$stability_theta$thetas))),
  ## censoring: 100 random panels per cohort
  num("random panels per cohort",            sprintf("%d random", cz$random_panels)),
  rule("q = 2 panel never above the 95th percentile", cz$cohorts_above_random_p95$survfcbf_q2_K10 == 0 &&
         cz$metabric_percentile$survfcbf_q2_K10 <= 0.95, sprintf("0 cohorts, METABRIC %.2f", cz$metabric_percentile$survfcbf_q2_K10)),
  wrd("finer binning above the 95th",        max(cz$cohorts_above_random_p95$survfcbf_q4_K10, cz$cohorts_above_random_p95$survfcbf_q10_K10),
      "does in %s of eleven"),
  ## simulation at matched size (sim_matched_size.R)
  rule("simulated sizes are three and ten",  identical(as.integer(s$sim_matched$K), c(3L, 10L)), "K = 3, 10"),
  num("sim recovery K=3, top-K",             fmt(sm2$survfcbf[["3"]], 2)),
  num("sim recovery K=10, top-K",            fmt(sm2$survfcbf[["10"]], 2)),
  num("sim recovery K=3, best default rival", fmt(sm2$best_rival[["3"]], 2)),
  num("sim recovery K=10, best default rival", fmt(sm2$best_rival[["10"]], 2)),
  num("sim recovery K=3, mRMR q=2",          fmt(sm2$mrmr_q2[["3"]], 2)),
  num("sim recovery K=10, mRMR q=2",         fmt(sm2$mrmr_q2[["10"]], 2)),
  rule("mRMR q=2 almost no copies",          all(unlist(s$sim_matched$redundant_copies[["3"]]$mrmr_q2,
                                                       s$sim_matched$redundant_copies[["10"]]$mrmr_q2) < 0.05),
       sprintf("%.3f / %.3f", s$sim_matched$redundant_copies[["3"]]$mrmr_q2, s$sim_matched$redundant_copies[["10"]]$mrmr_q2)),
  rule("top-K never selects two of a cluster", s$sim_matched$survfcbf_redundant_copies_total == 0, "0 copies"),
  ## gene-set enrichment moved to Supplementary Section S6 (its numbers are macros there)
  ## quoted twice, so each copy is checked in its phrase
  num("identical selections, abstract",      sprintf("selections on %d", s$validation$selections_identical)),
  num("identical selections, of generated",  sprintf("%d of %d", s$validation$selections_identical,
                                                     s$validation$selections_identical + s$validation$selections_skipped_exact_ties)),
  num("mismatches in the sweep with ties",   sprintf("all %d mismatches", s$validation$tie_sweep_mismatches)),
  num("datasets in the sweep with ties",     sprintf("sweep of %d datasets", s$validation$tie_sweep_datasets)),
  rule("no mismatch away from exact ties",   s$validation$selections_mismatches == 0 &&
         s$validation$tie_sweep_mismatches == s$validation$tie_sweep_mismatches_at_exact_tie,
       sprintf("%d of %d at a tie", s$validation$tie_sweep_mismatches_at_exact_tie, s$validation$tie_sweep_mismatches)),
  rule("SU matches the reference to 1e-15",  s$validation$su_max_abs_diff < 5e-15,
       format(s$validation$su_max_abs_diff, digits = 3)),
  num("single-cell dataset genes",           big(s$validation$scdengue_genes)),
  wrd("single-cell thresholds",              s$validation$scdengue_thresholds, "%s thresholds"),
  rule("single-cell selections identical",   s$validation$scdengue_mismatches == 0, "0 mismatches"),
  num("top-K exactness comparisons",         big(s$validation$topk_comparisons)),
  rule("top-K exact in every comparison",    s$validation$topk_mismatches == 0, "0 mismatches"),
  ## sourced since benchmark_metabric.csv (lambda.min, sha256 46bee94b...) is tracked and
  ## make_statistics.R reads it into stats$metabric
  num("METABRIC features",                   big(s$metabric$p)),
  num("METABRIC patients",                   big(s$metabric$n)),
  num("METABRIC ranker panel size",          big(s$metabric$ranker_panel)))


check_doc <- function(label, path, checks) {
  tex <- read_tex(path)
  cat(sprintf("=== %s number check (%s) ===\n", label, paste(basename(path), collapse = ", ")))
  ok <- vapply(checks, found, logical(1), tex = tex)
  for (i in seq_along(checks))
    cat(sprintf("  %-36s %-22s %s\n", checks[[i]]$name, checks[[i]]$value,
                if (ok[i]) "ok" else if (checks[[i]]$type == "rule") "FAILS" else "MISSING"))
  cat(sprintf("%d of %d claims found\n\n", sum(ok), length(ok)))
  sum(!ok)
}

bad <- check_doc("Application Note",
                 file.path("paper/appnote", c("appnote.tex", "abstract.tex", "body.tex")),
                 appnote_checks)

if (bad) {
  cat(sprintf("%d claim(s) absent or contradicted. Values absent from the text are either not\n", bad))
  cat("yet written or no longer match the results. Check each one before treating the\n")
  cat("document as current.\n")
  quit(status = 1L)
}
cat("All quoted figures trace to paper/results/paper_statistics.json\n")
