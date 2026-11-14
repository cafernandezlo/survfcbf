# Tables and quoted figures of the supplementary material, generated from the tracked
# result files so that no number is copied by hand. Each table is a LaTeX fragment in
# paper/appnote/generated/, read by paper/appnote/supplement_body.tex, and every figure
# quoted in the supplementary prose is a macro in generated/supp_numbers.tex.
#
# The values derived here are also checked against paper/results/paper_statistics.json,
# the source of the numbers in the Application Note, so the note and its supplementary
# material cannot drift apart. Run paper/make_statistics.R first.
#
# Usage (from the repository root):
#   Rscript --no-init-file paper/make_supplementary_tables.R
suppressWarnings(suppressMessages(library(jsonlite)))
root <- Sys.getenv("SURVFCBF_ROOT", unset = ".")
setwd(root)
RES <- "experiments/benchmark/results"
OUT <- "paper/appnote/generated"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
js <- jsonlite::fromJSON("paper/results/paper_statistics.json")

## ---- formatting -----------------------------------------------------------------------
## fixed decimals with a LaTeX thousands separator, "0.000" rather than "-0.000"
num <- function(x, d = 3) {
  s <- formatC(as.numeric(x), format = "f", digits = d, big.mark = ",")
  s <- gsub(",", "{,}", s, fixed = TRUE)
  s <- sub("^-(0\\.0+)$", "\\1", s)
  s <- sub("^-", "$-$", s)                       # a minus sign, not a hyphen
  s[is.na(x)] <- NA_character_
  s
}
int <- function(x) num(round(as.numeric(x)), 0)
pct <- function(x, d = 1) num(100 * as.numeric(x), d)
tr  <- function(...) paste0(paste(c(...), collapse = " & "), " \\\\")
word_of <- function(n) c("none", "one", "two", "three", "four", "five", "six", "seven", "eight",
                         "nine", "ten", "eleven")[as.integer(n) + 1L]
put <- function(name, lines) writeLines(lines, file.path(OUT, paste0(name, ".tex")))
## every value shared with the note must round to the JSON value
agree <- function(what, mine, theirs, d = 3) {
  if (!isTRUE(all(round(as.numeric(mine), d) == round(as.numeric(theirs), d))))
    stop(sprintf("%s: supplementary gives %s, paper_statistics.json %s", what,
                 paste(round(mine, d), collapse = ", "), paste(theirs, collapse = ", ")),
         call. = FALSE)
}

## The tuned-q variant of the real-data benchmark is left out: it was not rerun with the
## Spearman tie-break, so its rows would not reproduce with the package (Section S10). Its
## rows stay in benchmark_realdata.csv.
SURV <- c(survfcbf_topk = "top-$K$", survfcbf = "permutation")
BASE <- c(variance = "Variance", carss = "CAR", mim_surv = "MIM", cox_glmnet_l1 = "Lasso Cox")
SEL  <- c(SURV, BASE)
## seven columns of "0.576 (0.018)" need narrower gaps to fit the text width
head_sel <- function(first = "Cohort") c(
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{@{}l rr rrrr@{}}", "\\toprule",
  tr("", "\\multicolumn{2}{c}{\\pkg{survfcbf}}", "\\multicolumn{4}{c}{Baselines}"),
  "\\cmidrule(lr){2-3}\\cmidrule(l){4-7}",
  tr(first, SEL), "\\midrule")
foot <- c("\\bottomrule", "\\end{tabular}")

TYPE <- c(BLCA = "Bladder urothelial carcinoma", BRCA = "Breast invasive carcinoma",
          HNSC = "Head and neck squamous cell carcinoma", KIRC = "Kidney renal clear cell carcinoma",
          LGG = "Brain lower grade glioma", LUAD = "Lung adenocarcinoma",
          LUSC = "Lung squamous cell carcinoma", OV = "Ovarian serous cystadenocarcinoma",
          PAAD = "Pancreatic adenocarcinoma", SKCM = "Skin cutaneous melanoma",
          STAD = "Stomach adenocarcinoma")

## ---- inputs ---------------------------------------------------------------------------
rd <- utils::read.csv(file.path(RES, "benchmark_realdata.csv"), check.names = FALSE)
mb <- utils::read.csv(file.path(RES, "realdata_metabric/benchmark_metabric.csv"), check.names = FALSE)
COH <- sort(unique(rd$dataset))
stopifnot("unexpected cohorts" = setequal(COH, names(TYPE)))
cf <- list.files(file.path(RES, "realdata_cardinality"), pattern = "^[A-Z]+\\.csv$", full.names = TRUE)
cd <- do.call(rbind, lapply(cf, utils::read.csv, check.names = FALSE))
bommert <- new.env()
sys.source("experiments/benchmark/realdata/load_bommert.R", envir = bommert)
grid <- new.env()
sys.source("experiments/benchmark/grid_lib.R", envir = grid)
macros <- list()
mac <- function(name, value) macros[[name]] <<- value

## ---- Table: cohorts -------------------------------------------------------------------
coh <- unique(rd[, c("dataset", "n", "events", "p")])
coh <- coh[match(COH, coh$dataset), ]
pcod <- unique(cd[, c("cohort", "p")])
stopifnot("cardinality cohorts differ" = setequal(pcod$cohort, COH))
agree("TCGA patients", sum(coh$n), js$cohorts$total_patients, 0)
agree("TCGA events", sum(coh$events), js$cohorts$total_events, 0)
mbc <- unique(mb[, c("n", "events", "p")])
stopifnot(nrow(mbc) == 1L)
agree("METABRIC n, events, p", unlist(mbc), c(js$metabric$n, js$metabric$events, js$metabric$p), 0)
lines <- c("\\begin{tabularx}{\\linewidth}{@{}l >{\\raggedright\\arraybackslash}X l rrrrr@{}}", "\\toprule",
           tr("Cohort", "Cancer type", "Data source", "$n$", "Events", "Censored (\\%)", "$p$",
              "$p$ coding"), "\\midrule")
for (i in seq_len(nrow(coh))) {
  co <- coh$dataset[i]
  lines <- c(lines, tr(co, TYPE[[co]], paste(bommert$BOMMERT_IDS[[co]], collapse = ", "),
                       int(coh$n[i]), int(coh$events[i]), pct(1 - coh$events[i] / coh$n[i]),
                       int(coh$p[i]), int(pcod$p[pcod$cohort == co])))
}
lines <- c(lines, "\\midrule",
           tr("TCGA", "Eleven cohorts", "", int(sum(coh$n)), int(sum(coh$events)),
              pct(1 - sum(coh$events) / sum(coh$n)), "", ""),
           tr("METABRIC", "Breast cancer", "cBioPortal", int(mbc$n), int(mbc$events),
              pct(1 - mbc$events / mbc$n), int(mbc$p), "n.r."),
           "\\bottomrule", "\\end{tabularx}")
put("tab_cohorts", lines)
mac("SuppTcgaPatients", int(sum(coh$n))); mac("SuppTcgaEvents", int(sum(coh$events)))
mac("SuppCodingMin", int(min(pcod$p))); mac("SuppCodingMax", int(max(pcod$p)))
mac("SuppMetabricN", int(mbc$n)); mac("SuppMetabricEvents", int(mbc$events)); mac("SuppMetabricP", int(mbc$p))

## ---- Tables: per-cohort results -------------------------------------------------------
cell <- function(d, co, s, metric, digits, se, star) {
  z <- d[d$dataset == co & d$selector == s, , drop = FALSE]
  if (!nrow(z)) return("n.r.")
  v <- z[[metric]]
  if (is.na(v)) return("n.m.")
  out <- num(v, digits)
  if (se && !is.na(z[[paste0(metric, "_se")]])) out <- sprintf("%s (%s)", out, num(z[[paste0(metric, "_se")]], digits))
  if (star && z$n_failed > 0) out <- paste0(out, "$^{\\ast}$")
  out
}
per_cohort <- function(name, metric, digits = 3, se = TRUE, star = TRUE, json = NULL) {
  lines <- head_sel()
  for (co in COH) lines <- c(lines, tr(co, vapply(names(SEL), function(s) cell(rd, co, s, metric, digits, se, star), "")))
  m <- vapply(names(SEL), function(s) mean(rd[[metric]][rd$selector == s], na.rm = TRUE), 0)
  k <- vapply(names(SEL), function(s) sum(!is.na(rd[[metric]][rd$selector == s])), 0)
  if (!is.null(json)) agree(paste("mean", metric), m, unlist(json[names(SEL)]), digits)
  lines <- c(lines, "\\midrule", tr("Mean", ifelse(is.nan(m), "n.m.", num(m, digits))),
             tr("Cohorts", as.character(k)), "\\midrule",
             tr("METABRIC", vapply(names(SEL), function(s) cell(mb, "METABRIC", s, metric, digits, se, star), "")),
             foot)
  put(name, lines)
}
per_cohort("tab_unoc", "uno_c", json = js$benchmark$uno_c)
per_cohort("tab_ibs", "ibs")
per_cohort("tab_stability", "stability", se = FALSE, star = FALSE, json = js$benchmark$stability)
per_cohort("tab_panel", "parsimony", digits = 1, se = FALSE, star = FALSE, json = js$benchmark$parsimony)

lines <- head_sel()
fail_cell <- function(d, co, s) {
  z <- d[d$dataset == co & d$selector == s, , drop = FALSE]
  if (!nrow(z)) "n.r." else sprintf("%d (%d)", as.integer(z$n_failed), as.integer(z$n_null_fit))
}
for (co in COH) lines <- c(lines, tr(co, vapply(names(SEL), function(s) fail_cell(rd, co, s), "")))
tot <- vapply(names(SEL), function(s) sprintf("%d (%d)", sum(rd$n_failed[rd$selector == s]),
                                            sum(rd$n_null_fit[rd$selector == s])), "")
lines <- c(lines, "\\midrule", tr("Total", tot), "\\midrule",
           tr("METABRIC", vapply(names(SEL), function(s) fail_cell(mb, "METABRIC", s), "")), foot)
put("tab_failures", lines)
agree("null-fit folds", vapply(names(SEL), function(s) sum(rd$n_null_fit[rd$selector == s]), 0),
      unlist(js$benchmark$null_fit_folds[names(SEL)]), 0)

## ---- Table: Bayesian signed-rank against top-K ----------------------------------------
bsr <- utils::read.csv(file.path(RES, "realdata_crossdataset/bayes_signed_rank.csv"), check.names = FALSE)
RIV <- c("variance", "carss", "mim_surv", "cox_glmnet_l1", "survfcbf")
stopifnot(all(RIV %in% bsr$rival))
bsr <- bsr[bsr$rival %in% RIV, ]
mat <- function(metric) tapply(rd[[metric]], list(rd$dataset, rd$selector), identity)
Ms <- mat("stability")
pp <- function(p) ifelse(p > 0.95, sprintf("\\textbf{%s}", num(p, 2)), num(p, 2))
lines <- c("\\begin{tabular}{@{}l rrrr rrrr@{}}", "\\toprule",
           tr("", "\\multicolumn{4}{c}{Uno's C, ROPE $\\pm 0.02$}",
              "\\multicolumn{4}{c}{Stability, ROPE $\\pm 0.05$}"),
           "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}",
           tr("Rival", "Cohorts", "top-$K$", "ROPE", "Rival", "Cohorts", "top-$K$", "ROPE", "Rival"),
           "\\midrule")
for (rv in RIV) {
  z <- bsr[bsr$rival == rv, ]
  nS <- sum(is.finite(Ms[, "survfcbf_topk"]) & is.finite(Ms[, rv]))
  lab <- if (rv %in% names(BASE)) BASE[[rv]] else paste("\\pkg{survfcbf},", SURV[[rv]])
  lines <- c(lines, tr(lab, z$C_n, pp(z$C_topk), pp(z$C_rope), pp(z$C_rival),
                       nS, pp(z$S_topk), pp(z$S_rope), pp(z$S_rival)))
}
put("tab_bayes", c(lines, foot))
## ---- Repeated cross-validation, the primary analysis (S4) ------------------------------
## Cohort means over ten five-fold splits of the truncated Uno's C (realdata/repeated_cv_main.R
## for the headline, realdata/repeated_cv.R for the matched sizes)
RK <- c("variance", "carss", "mim_surv")
rcv_read <- function(dir) do.call(rbind, lapply(list.files(file.path(RES, dir), pattern = "^[A-Z]+_rep[0-9]+\\.csv$",
                                                           full.names = TRUE), utils::read.csv, check.names = FALSE))
rcv_inf <- function(dir) { z <- utils::read.csv(file.path(RES, dir, "inference.csv"), check.names = FALSE)
  z$rival <- sub("_K[0-9]+$", "", z$rival); z }
rmn <- rcv_read("realdata_repeated_cv_main")
cmn <- aggregate(cbind(uno_c_trunc, uno_c, stability) ~ cohort + selector, rmn, mean)
Mt <- tapply(cmn$uno_c_trunc, list(cmn$cohort, cmn$selector), identity)
## 95% t intervals of the paired gap, top-K minus each ranker, over the cohort means
gap_ci <- lapply(stats::setNames(nm = RK), function(r) { d <- Mt[, "survfcbf_topk"] - Mt[, r]
  c(mean(d), stats::t.test(d)$conf.int) })
agree("gap CIs vs JSON", unlist(gap_ci), unlist(lapply(RK, function(r) unlist(js$repeated_cv$headline$gap_ci_trunc[[r]]))))
for (r in RK) {
  nm <- c(variance = "Var", carss = "Car", mim_surv = "Mim")[[r]]
  mac(paste0("SuppGap", nm), num(gap_ci[[r]][1])); mac(paste0("SuppGapLo", nm), num(gap_ci[[r]][2]))
  mac(paste0("SuppGapHi", nm), num(gap_ci[[r]][3]))
}
## the same comparisons under the mean-mass summary of the posterior
rsh <- utils::read.csv(file.path(RES, "realdata_repeated_cv_main/rope_sensitivity.csv"), check.names = FALSE)
pmean <- function(r) rsh$p_rival[rsh$rival == r & rsh$metric == "uno_c_trunc" & rsh$rope == 0.02 & rsh$summary == "mean"]
mac("SuppMeanVar", num(pmean("variance"), 2)); mac("SuppMeanCar", num(pmean("carss"), 2))
mac("SuppMeanMim", num(pmean("mim_surv"), 2))
n_dec <- sum(vapply(RK, pmean, 0) > 0.95)
mac("SuppMeanDecidedText", c("none of the three comparisons would be decided",
                             "only one of the three comparisons would be decided",
                             "two of the three comparisons would be decided",
                             "all three comparisons would still be decided")[n_dec + 1L])
hinf <- rcv_inf("realdata_repeated_cv_main")
VSHORT <- c("survfcbf better" = "\\pkg{survfcbf}", "rival better" = "rival",
            "practically equivalent" = "equivalent", "undecided at n=11" = "undecided")
MET3 <- c(uno_c_trunc = "Uno's C, truncated", uno_c = "Uno's C, untruncated", stability = "Stability")
lines <- c("\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}l l rrr r rrr r@{}}", "\\toprule",
           tr("Measure", "Rival", "\\pkg{survfcbf}", "Rival", "Difference (SE)", "Cohorts", "\\pkg{survfcbf}",
              "ROPE", "Rival", "Splits"), "\\midrule")
for (m in names(MET3)) {
  for (r in RK) {
    z <- hinf[hinf$rival == r & hinf$metric == m, ]
    lines <- c(lines, tr(if (r == RK[1]) MET3[[m]] else "", BASE[[r]], num(z$focus_mean), num(z$rival_mean),
                         sprintf("%s (%s)", num(z$mean_diff), num(z$se_diff)),
                         sprintf("%d/%d", z$wins_focus, z$wins_rival), pp(z$p_focus), pp(z$p_rope), pp(z$p_rival),
                         sprintf("%d/%d", z$reps_decided_same, z$reps)))
  }
  lines <- c(lines, "\\addlinespace")
}
## integrated Brier score, lower is better, without a ROPE and therefore without the test
cmi <- aggregate(ibs ~ cohort + selector, rmn, mean)
Mi <- tapply(cmi$ibs, list(cmi$cohort, cmi$selector), identity)
ibs_gap <- list()
for (r in RK) {
  d <- Mi[, "survfcbf_topk"] - Mi[, r]; ibs_gap[[r]] <- c(mean(d), sum(d > 0))
  agree(paste("IBS gap", r), mean(d), js$repeated_cv$headline$gap_ci_ibs[[r]]$diff)
  lines <- c(lines, tr(if (r == RK[1]) "Integrated Brier score" else "", BASE[[r]], num(mean(Mi[, "survfcbf_topk"])),
                       num(mean(Mi[, r])), sprintf("%s (%s)", num(mean(d)), num(stats::sd(d) / sqrt(length(d)))),
                       sprintf("%d/%d", sum(d < 0), sum(d > 0)), "", "", "", ""))
}
put("tab_rcv", c(lines, foot))
ig <- vapply(ibs_gap, function(v) v[1], 0); iw <- vapply(ibs_gap, function(v) v[2], 0)
stopifnot("every IBS interval includes zero, as the text says" =
            all(vapply(RK, function(r) js$repeated_cv$headline$gap_ci_ibs[[r]]$lo < 0 &&
                         js$repeated_cv$headline$gap_ci_ibs[[r]]$hi > 0, TRUE)),
          "top-K has the higher IBS against every ranker" = all(ig > 0))
mac("SuppIbsGapMin", num(min(ig))); mac("SuppIbsGapMax", num(max(ig)))
mac("SuppIbsWorseMin", word_of(min(iw))); mac("SuppIbsWorseMax", word_of(max(iw)))
mac("SuppRcvSplits", word_of(js$repeated_cv$n_splits))
## per-cohort truncated C, ten-split means, and METABRIC
mrc <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_metabric/repeated_cv"), pattern = "^rep[0-9]+\\.csv$",
                                        full.names = TRUE), utils::read.csv, check.names = FALSE))
HS <- c(survfcbf_topk = "top-$K$", variance = "Variance", carss = "CAR", mim_surv = "MIM")
csd <- aggregate(uno_c_trunc ~ cohort + selector, rmn, stats::sd)
lines <- c("\\begin{tabular}{@{}l rrrr@{}}", "\\toprule", tr("Cohort", HS), "\\midrule")
cellr <- function(co, s, d, sdd) sprintf("%s (%s)", num(d$uno_c_trunc[d$cohort == co & d$selector == s]),
                                        num(sdd$uno_c_trunc[sdd$cohort == co & sdd$selector == s]))
for (co in COH) lines <- c(lines, tr(co, vapply(names(HS), function(s) cellr(co, s, cmn, csd), "")))
mm <- vapply(names(HS), function(s) mean(cmn$uno_c_trunc[cmn$selector == s]), 0)
agree("headline truncated means", mm, unlist(js$repeated_cv$headline$means$uno_c_trunc[names(HS)]))
mrm <- aggregate(uno_c_trunc ~ selector, mrc, mean); mrs <- aggregate(uno_c_trunc ~ selector, mrc, stats::sd)
mrm$cohort <- mrs$cohort <- "METABRIC"
lines <- c(lines, "\\midrule", tr("Mean", num(mm)), "\\midrule",
           tr("METABRIC", vapply(names(HS), function(s) cellr("METABRIC", s, mrm, mrs), "")), foot)
put("tab_rcv_cohort", lines)
mg <- js$repeated_cv$metabric$gap_trunc
mac("SuppMetabricGapMin", num(min(-vapply(mg, function(g) g$diff, 0))))
mac("SuppMetabricGapMax", num(max(-vapply(mg, function(g) g$diff, 0))))
mac("SuppMetabricGlmnetDiff", num(js$repeated_cv$metabric$split1_max_abs_c_diff_vs_published))
mgi <- js$repeated_cv$metabric$gap_ibs
stopifnot("METABRIC IBS worse for top-K in every split" = all(vapply(mgi, function(g) g$splits_topk_worse == 10L, TRUE)))
mac("SuppMetabricIbsMin", num(min(vapply(mgi, function(g) g$diff, 0))))
mac("SuppMetabricIbsMax", num(max(vapply(mgi, function(g) g$diff, 0))))
## clinical reference on METABRIC (realdata/metabric_clinical.R), published split
cl <- utils::read.csv(file.path(RES, "realdata_metabric/clinical_reference.csv"), check.names = FALSE)
clv <- function(m, col) cl[[col]][cl$model == m]
agree("clinical increment", clv("clinical_plus_panel_minus_clinical", "uno_c_trunc"), js$metabric_clinical$added_by_panel_trunc)
mac("SuppClinN", int(clv("clinical", "n"))); mac("SuppClinEvents", int(clv("clinical", "events")))
mac("SuppClinC", num(clv("clinical", "uno_c_trunc"))); mac("SuppClinPanelC", num(clv("panel", "uno_c_trunc")))
mac("SuppClinBothC", num(clv("clinical_plus_panel", "uno_c_trunc")))
mac("SuppClinAdd", num(clv("clinical_plus_panel_minus_clinical", "uno_c_trunc")))
mac("SuppClinAddSe", num(clv("clinical_plus_panel_minus_clinical", "uno_c_trunc_se")))
## sensitivity to the ROPE and to the posterior summary: headline and matched K = 10
rs10 <- utils::read.csv(file.path(RES, "realdata_repeated_cv/rope_sensitivity.csv"), check.names = FALSE)
rs10$rival <- sub("_K[0-9]+$", "", rs10$rival)
RV6 <- c(variance = "Variance", carss = "CAR", mim_surv = "MIM", mrmr = "mRMR",
         mim_surv_q2 = "MIM, $q = 2$", mrmr_q2 = "mRMR, $q = 2$")
vcell <- function(z, r, m, rope, summ = "argmax") VSHORT[[z$verdict[z$rival == r & z$metric == m & z$rope == rope & z$summary == summ]]]
rope_row <- function(z, r, lab) tr(lab, vcell(z, r, "uno_c_trunc", 0.01), vcell(z, r, "uno_c_trunc", 0.02),
                                   vcell(z, r, "uno_c_trunc", 0.03), vcell(z, r, "uno_c_trunc", 0.02, "mean"),
                                   vcell(z, r, "stability", 0.03), vcell(z, r, "stability", 0.05),
                                   vcell(z, r, "stability", 0.07))
lines <- c("\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}l llll lll@{}}", "\\toprule",
           tr("", "\\multicolumn{4}{c}{Uno's C, truncated}", "\\multicolumn{3}{c}{Stability}"),
           "\\cmidrule(lr){2-5}\\cmidrule(l){6-8}",
           tr("Rival", "$\\pm 0.01$", "$\\pm 0.02$", "$\\pm 0.03$", "$\\pm 0.02$, mean", "$\\pm 0.03$", "$\\pm 0.05$", "$\\pm 0.07$"),
           "\\midrule", "\\multicolumn{8}{@{}l}{\\textit{A. Headline, at the rankers' own panel size}} \\\\")
for (r in RK) lines <- c(lines, rope_row(rsh, r, BASE[[r]]))
lines <- c(lines, "\\midrule", "\\multicolumn{8}{@{}l}{\\textit{B. Matched panel size, $K = 10$}} \\\\")
for (r in names(RV6)) lines <- c(lines, rope_row(rs10, r, RV6[[r]]))
put("tab_rope", c(lines, foot))

## ---- Table: power ---------------------------------------------------------------------
## two families of comparisons (realdata/power_analysis.R): the headline (top-K against the
## three rankers at their cut-off) and the matched sizes (the same rankers at K = 5, 10, 20),
## both on the ten-split cohort means of the truncated C
pw <- utils::read.csv(file.path(RES, "realdata_power/pairwise_sd.csv"), check.names = FALSE)
sg2 <- utils::read.csv(file.path(RES, "realdata_power/sigma.csv"), check.names = FALSE)
an <- utils::read.csv(file.path(RES, "realdata_power/analytic_n.csv"), check.names = FALSE)
sp <- utils::read.csv(file.path(RES, "realdata_power/simulated_power.csv"), check.names = FALSE)
REPS <- unique(sp$reps); stopifnot("one replicate count" = length(REPS) == 1L)
FAM <- c(headline = "Headline", matched = "Matched size")
sgf <- function(f, col) sg2[[col]][sg2$family == f]
for (f in names(FAM)) {
  agree(paste("sigma", f), sgf(f, "sigma"), js$power$sigma[[f]]$sigma)
  agree(paste("sigma range", f), c(sgf(f, "sd_min"), sgf(f, "sd_max")),
        c(js$power$sigma[[f]]$sd_min, js$power$sigma[[f]]$sd_max))
  stopifnot("pairs per family" = sgf(f, "n_pairs") == sum(pw$family == f))
}
nm_f <- c(headline = "Head", matched = "Matched")
for (f in names(FAM)) {
  mac(paste0("SuppSigma", nm_f[[f]]), num(sgf(f, "sigma")))
  mac(paste0("SuppSigma", nm_f[[f]], "Min"), num(sgf(f, "sd_min")))
  mac(paste0("SuppSigma", nm_f[[f]], "Max"), num(sgf(f, "sd_max")))
  mac(paste0("SuppPairs", nm_f[[f]]), word_of(sgf(f, "n_pairs")))
  mac(paste0("SuppDelta", nm_f[[f]]), num(sgf(f, "delta_80_n11")))
  mac(paste0("SuppNeedThree", nm_f[[f]]), int(an$n_80[an$family == f & an$delta == 0.03]))
}
mac("SuppNeedThreeMatchedMax", int(an$n_80_pessimistic[an$family == "matched" & an$delta == 0.03]))
mx <- pw[pw$family == "matched", ][which.max(pw$sd_diff[pw$family == "matched"]), ]
mac("SuppSigmaMatchedMaxWhat", sprintf("against the %s filter at $K = %s$",
                                       if (grepl("^variance", mx$b) || grepl("^variance", mx$a)) "variance" else "CAR",
                                       sub(".*_K", "", mx$a)))
mac("SuppSdLegacy", num(stats::median(pw$sd_diff[pw$family == "legacy"])))
agree("legacy SD", stats::median(pw$sd_diff[pw$family == "legacy"]), js$power$legacy_sd_median)
mac("SuppPowerReps", int(REPS))
lines <- c("\\begin{tabular}{@{}l rrr rrr@{}}", "\\toprule",
           "\\multicolumn{7}{@{}l}{\\textit{A. Cohorts needed by a paired test (two-sided $\\alpha = 0.05$)}} \\\\",
           tr("", sprintf("\\multicolumn{3}{c}{%s, SD %s}", FAM[["headline"]], num(sgf("headline", "sigma"))),
              sprintf("\\multicolumn{3}{c}{%s, SD %s}", FAM[["matched"]], num(sgf("matched", "sigma")))),
           "\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
           tr("Difference", rep(c("80\\%", "90\\%", "80\\%, largest SD"), 2)), "\\midrule")
for (dl in sort(unique(an$delta))) {
  a <- function(f, col) int(an[[col]][an$family == f & an$delta == dl])
  lines <- c(lines, tr(num(dl, 2), unlist(lapply(names(FAM), function(f)
    c(a(f, "n_80"), a(f, "n_90"), a(f, "n_80_pessimistic"))))))
}
Ns <- sort(unique(sp$N))
lines <- c(lines, "\\midrule",
           "\\multicolumn{7}{@{}l}{\\textit{B. Power of the 0.95 rule (Bayesian signed-rank, ROPE $\\pm 0.02$)}} \\\\",
           tr("Difference", rep(sprintf("%d cohorts", Ns), 2)), "\\midrule")
for (dl in sort(unique(sp$delta))) {
  cellp <- function(f) vapply(Ns, function(N) { r <- sp$rate_p95[sp$family == f & sp$N == N & sp$delta == dl]
    sprintf("%s (%s)", num(r, 2), num(sqrt(r * (1 - r) / REPS), 2)) }, "")
  lines <- c(lines, tr(num(dl, 2), cellp("headline"), cellp("matched")))
}
put("tab_power", c(lines, foot))
mac("SuppPowerHeadThree", num(sp$rate_p95[sp$family == "headline" & sp$N == 11 & sp$delta == 0.03], 2))
mac("SuppPowerHeadFive", num(sp$rate_p95[sp$family == "headline" & sp$N == 11 & sp$delta == 0.05], 2))

## ---- Tables: matched panel size -------------------------------------------------------
CARD <- c(survfcbf = "\\pkg{survfcbf}", variance = "Variance", carss = "CAR", mim_surv = "MIM")
Ks <- sort(unique(cd$K))
panel <- function(title, f) {
  out <- c(sprintf("\\multicolumn{5}{@{}l}{\\textit{%s}} \\\\", title))
  for (K in Ks) out <- c(out, tr(int(K), vapply(names(CARD), function(b) {
    z <- cd[cd$base == b & cd$K == K, ]
    if (!nrow(z)) "n.r." else f(z) }, "")))
  out
}
mean_c <- function(z) num(mean(z$uno_c, na.rm = TRUE))
mean_s <- function(z) num(mean(z$stability, na.rm = TRUE))
for (b in names(CARD)) for (K in Ks[Ks <= max(cd$K[cd$base == b])]) {
  z <- cd[cd$base == b & cd$K == K, ]
  agree(sprintf("C at K=%d, %s", K, b), mean(z$uno_c, na.rm = TRUE), js$cardinality$uno_c_by_K[[b]][[paste0("K", K)]])
  agree(sprintf("stability at K=%d, %s", K, b), mean(z$stability, na.rm = TRUE),
        js$cardinality$stability_by_K[[b]][[paste0("K", K)]])
}
stopifnot("matched-size runs lost folds" = all(cd$n_failed == 0))
## top-K returns fewer than K features when the whole feature set holds fewer than K
## predominant features, so the survfcbf panel is only matched to K up to a point
size_rng <- function(z) if (all(abs(z$parsimony - z$K) < 1e-9)) int(z$K[1]) else
  sprintf("%s (%s--%s)", num(mean(z$parsimony), 1), num(min(z$parsimony), 1), num(max(z$parsimony), 1))
fc <- cd[cd$base == "survfcbf", ]
stopifnot("rankers must sit exactly at K" = all(abs(cd$parsimony[cd$base != "survfcbf"] -
                                                    cd$K[cd$base != "survfcbf"]) < 1e-9))
short <- vapply(sort(unique(fc$K)), function(K) sum(fc$parsimony[fc$K == K] < K - 1e-9), 0)
names(short) <- sort(unique(fc$K))
mac("SuppCardExactMax", int(max(as.integer(names(short))[short == 0])))
mac("SuppCardShortTwenty", int(short[["20"]]))
mac("SuppCardPanelMin", num(min(fc$parsimony[fc$K == max(fc$K)]), 1))
mac("SuppCardPanelMax", num(max(fc$parsimony[fc$K == max(fc$K)]), 1))
mac("SuppCardPanelMean", num(mean(fc$parsimony[fc$K == max(fc$K)]), 1))
## matched sizes over ten splits (realdata/repeated_cv.R at K = 5, 10, 20)
MKD <- c("5" = "realdata_repeated_cv_K5", "10" = "realdata_repeated_cv", "20" = "realdata_repeated_cv_K20")
minf <- do.call(rbind, lapply(names(MKD), function(K) cbind(K = as.integer(K), rcv_inf(MKD[[K]]))))
rng_m <- function(r) range(minf$mean_diff[minf$rival == r & minf$metric == "uno_c_trunc"])
for (r in names(RV6)) agree(paste("matched range", r), rng_m(r), js$repeated_cv$matched$range_trunc[[r]])
mac("SuppRcvVarLo", num(rng_m("variance")[1])); mac("SuppRcvVarHi", num(rng_m("variance")[2]))
mac("SuppRcvCarLo", num(rng_m("carss")[1])); mac("SuppRcvCarHi", num(rng_m("carss")[2]))
mac("SuppRcvMimLo", num(-rng_m("mim_surv")[2])); mac("SuppRcvMimHi", num(-rng_m("mim_surv")[1]))
nfs <- lapply(names(MKD), function(K) { r <- rcv_read(MKD[[K]]); r$sel <- sub("_K[0-9]+$", "", r$selector)
  data.frame(K = as.integer(K), sel = r$sel, null = r$n_null_fit) })
nfs <- do.call(rbind, nfs)
nshare <- tapply(nfs$null, list(nfs$K, nfs$sel), sum) / (5 * tapply(nfs$null, list(nfs$K, nfs$sel), length))
mac("SuppRcvNullVarLo", pct(min(nshare[, "variance"]), 0)); mac("SuppRcvNullVarHi", pct(max(nshare[, "variance"]), 0))
mac("SuppRcvNullTopk", int(sum(nfs$null[nfs$sel == "survfcbf_q2"])))
mac("SuppRcvNullTopkFolds", int(5L * sum(nfs$sel == "survfcbf_q2")))
lines <- c("\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}l r r r l r r l@{}}", "\\toprule",
           tr("", "", "\\multicolumn{3}{c}{Uno's C, truncated}", "\\multicolumn{3}{c}{Stability}"),
           "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
           tr("Rival", "$K$", "Difference (SE)", "Cohorts", "Verdict", "Difference (SE)", "Cohorts", "Verdict"),
           "\\midrule")
for (r in names(RV6)) {
  for (K in c(5L, 10L, 20L)) {
    zc <- minf[minf$rival == r & minf$K == K & minf$metric == "uno_c_trunc", ]
    zs <- minf[minf$rival == r & minf$K == K & minf$metric == "stability", ]
    lines <- c(lines, tr(if (K == 5L) RV6[[r]] else "", K,
                         sprintf("%s (%s)", num(zc$mean_diff), num(zc$se_diff)), sprintf("%d/%d", zc$wins_focus, zc$wins_rival),
                         VSHORT[[zc$verdict]],
                         sprintf("%s (%s)", num(zs$mean_diff), num(zs$se_diff)), sprintf("%d/%d", zs$wins_focus, zs$wins_rival),
                         VSHORT[[zs$verdict]]))
  }
  if (r != tail(names(RV6), 1)) lines <- c(lines, "\\addlinespace")
}
put("tab_rcv_matched", c(lines, foot))
## every 0.95-rule decision reported in the supplementary tables, for the multiplicity note
n_tab <- c(bayes = 2L * length(RIV), cardinality = nrow(utils::read.csv(file.path(RES,
             "realdata_cardinality/cardinality_inference.csv"))), rcv = 3L * length(RK),
           rcv_matched = 2L * 3L * length(RV6))
mac("SuppNComparisons", int(sum(n_tab)))
lines <- c("\\begin{tabular}{@{}r lrrr@{}}", "\\toprule", tr("$K$", CARD), "\\midrule",
           panel("A. Panel size, mean (range) over cohorts", size_rng), "\\midrule",
           panel("B. Uno's C, mean over cohorts", mean_c), "\\midrule",
           panel("C. Stability, mean over cohorts", mean_s), "\\midrule",
           panel(sprintf("D. Null-fit folds, out of %d", 5L * length(COH)),
                 function(z) as.character(sum(z$n_null_fit))), foot)
put("tab_cardinality", lines)

inf <- utils::read.csv(file.path(RES, "realdata_cardinality/cardinality_inference.csv"), check.names = FALSE)
lines <- c("\\begin{tabular}{@{}l r r r r rrr@{}}", "\\toprule",
           tr("Rival", "$K$", "\\pkg{survfcbf} size", "Difference (SE)", "Cohorts", "\\pkg{survfcbf}",
              "ROPE", "Rival"),
           "\\midrule")
for (rv in c("variance", "carss", "mim_surv")) {
  z <- inf[inf$rival == rv, ]; z <- z[order(z$K), ]
  for (i in seq_len(nrow(z)))
    lines <- c(lines, tr(if (i == 1L) CARD[[rv]] else "", z$K[i],
                         num(mean(fc$parsimony[fc$K == z$K[i]]), 1),
                         sprintf("%s (%s)", num(z$mean_diff[i]), num(z$se_diff[i])),
                         sprintf("%d/%d", z$wins_focus[i], z$wins_rival[i]),
                         pp(z$p_focus[i]), pp(z$p_rope[i]), pp(z$p_rival[i])))
  if (rv != "mim_surv") lines <- c(lines, "\\addlinespace")
}
stopifnot("a matched-size comparison is decided; the note says none is" =
            all(pmax(inf$p_focus, inf$p_rope, inf$p_rival) <= 0.95))
put("tab_cardinality_inference", c(lines, foot))

## ---- Table: gene-set enrichment -------------------------------------------------------
gs <- readRDS(file.path(RES, "realdata_gsea/gsea.rds"))
msig <- readRDS("experiments/benchmark/realdata/prior_knowledge/msigdb_H_C6.rds")
pcg <- readRDS("experiments/benchmark/realdata/prior_knowledge/protein_coding.rds")
GS <- c(survfcbf = "\\pkg{survfcbf}", variance = "Variance", carss = "CAR", mim_surv = "MIM")
sig <- function(co, s) sum(gs$padj[gs$cohort == co & gs$selector == s] < 0.05, na.rm = TRUE)
mp <- utils::read.csv(file.path(RES, "realdata_gsea/mapping.csv"), check.names = FALSE)
mp <- mp[match(COH, mp$cohort), ]
stopifnot("mapping table must cover every cohort" = identical(mp$cohort, COH),
          "mapping arithmetic" = all(mp$coding_features - mp$no_symbol - mp$duplicated_symbol == mp$ranked))
lines <- c("\\begin{tabular}{@{}l rrr r rrrr@{}}", "\\toprule",
           tr("", "\\multicolumn{3}{c}{Genes}", "", "\\multicolumn{4}{c}{Significant sets}"),
           "\\cmidrule(lr){2-4}\\cmidrule(l){6-9}",
           tr("Cohort", "Coding", "Lost", "Ranked", "Sets tested", GS), "\\midrule")
for (i in seq_along(COH)) {
  co <- COH[i]
  tested <- unique(vapply(names(GS), function(s) sum(gs$cohort == co & gs$selector == s), 0))
  lines <- c(lines, tr(co, int(mp$coding_features[i]), int(mp$no_symbol[i] + mp$duplicated_symbol[i]),
                       int(mp$ranked[i]), paste(tested, collapse = "/"),
                       vapply(names(GS), function(s) sig(co, s), 0)))
}
hits <- vapply(names(GS), function(s) sum(vapply(COH, function(co) sig(co, s) > 0, TRUE)), 0)
agree("GSEA cohorts", hits, unlist(js$enrichment$gsea_cohorts_enriched[names(GS)]), 0)
lines <- c(lines, "\\midrule", tr("Cohorts with a set", "", "", "", "", hits), foot)
put("tab_gsea", lines)
mac("SuppDistinctSuMin", int(min(mp$distinct_su))); mac("SuppDistinctSuMax", int(max(mp$distinct_su)))
mac("SuppNoSymbolMin", int(min(mp$no_symbol))); mac("SuppNoSymbolMax", int(max(mp$no_symbol)))
mac("SuppDupMin", int(min(mp$duplicated_symbol))); mac("SuppDupMax", int(max(mp$duplicated_symbol)))
mac("SuppGseaVariance", int(hits[["variance"]]))

## ---- Table: composition of the q = 2 target -------------------------------------------
tg <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_target"), "^[A-Z]+\\.csv$", full.names = TRUE),
                            utils::read.csv, check.names = FALSE))
tg <- tg[match(COH, tg$cohort), ]
stopifnot("target composition must cover every cohort" = identical(tg$cohort, COH),
          "target composition arithmetic" = all(tg$upper_events + tg$upper_censored == tg$upper_n),
          "the target must match the benchmark cohorts" = all(tg$n == coh$n & tg$events == coh$events))
## median follow-up of the cohort by reverse Kaplan-Meier (realdata/center_followup.R), days
cfu <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_center"), "^[A-Z]+\\.csv$", full.names = TRUE),
                             utils::read.csv, check.names = FALSE))
cfu <- cfu[match(COH, cfu$cohort), ]
stopifnot("centre table must match the benchmark cohorts" = identical(cfu$cohort, COH) && all(cfu$n == coh$n))
lines <- c("\\begin{tabular}{@{}l rrr rr rrr@{}}", "\\toprule",
           tr("", "", "", "\\multicolumn{3}{c}{Upper class}", "\\multicolumn{3}{c}{Follow-up (days)}"),
           "\\cmidrule(lr){4-6}\\cmidrule(l){7-9}",
           tr("Cohort", "$n$", "Events", "Size", "Censored", "Share (\\%)", "Upper", "Lower", "Cohort"), "\\midrule")
for (i in seq_len(nrow(tg)))
  lines <- c(lines, tr(tg$cohort[i], int(tg$n[i]), int(tg$events[i]), int(tg$upper_n[i]),
                       int(tg$upper_censored[i]), pct(tg$upper_censored[i] / tg$upper_n[i]),
                       int(tg$followup_censored_upper[i]), int(tg$followup_censored_lower[i]),
                       int(cfu$median_followup_reverse_km[i])))
put("tab_target", c(lines, foot))
sh <- tg$upper_censored / tg$upper_n
mac("SuppTargetShareMax", pct(max(sh), 0)); mac("SuppTargetShareMaxCohort", tg$cohort[which.max(sh)])
mac("SuppTargetHalfOrMore", int(sum(sh >= 0.5)))
mac("SuppMsigdbrVersion", msig$meta$msigdbr_version); mac("SuppMsigdbDate", msig$meta$date)
mac("SuppSetsHallmark", int(msig$meta$n_sets[["H"]])); mac("SuppSetsCsix", int(msig$meta$n_sets[["C6"]]))
mac("SuppSetsTotal", int(sum(msig$meta$n_sets)))
mac("SuppCodingIds", int(length(pcg$ensembl))); mac("SuppCodingDate", pcg$date)
mac("SuppGseaSurvfcbf", int(hits[["survfcbf"]]))

## ---- Table: simulation grid -----------------------------------------------------------
cfg <- grid$grid_config("cesga_core")
sg <- utils::read.csv(file.path(RES, "raw_cesga_core.csv"), check.names = FALSE)
## A_cens50 repeats A_n600 (same configuration, same seeds): see paper/make_statistics.R
sg <- sg[sg$scenario != "A_cens50", ]
conf <- function(sc) {
  v <- as.numeric(sub("^[A-Z]_[a-z]+", "", sc))
  list(n = if (grepl("^A_n", sc)) v else cfg$n_ref,
       cens = if (grepl("^A_cens", sc)) v / 100 else cfg$censor_ref,
       p = if (grepl("^B_p", sc)) v else cfg$p_ref)
}
SC <- c(sprintf("A_n%d", cfg$n_levels), sprintf("A_cens%02.0f", 100 * setdiff(cfg$censor_levels, cfg$censor_ref)),
        sprintf("B_p%d", cfg$p_levels))
stopifnot("simulation scenarios differ from grid_lib.R" = setequal(SC, unique(sg$scenario)))
share <- function(x) sprintf("%s (%s)", pct(mean(x)), pct(sqrt(mean(x) * (1 - mean(x)) / length(x))))
sim_row <- function(z, lab) {
  ev <- z$events[z$selector == "variance"]
  out <- lab
  for (s in c("survfcbf", "survfcbf_cvq")) {
    w <- z[z$selector == s, ]
    out <- c(out, share(is.na(w$sensitivity)), share(w$n_failed >= cfg$outer))
  }
  c(out[1:3], sprintf("%s (%s)", num(mean(ev), 0), num(stats::sd(ev) / sqrt(length(ev)), 1)), out[-(1:3)])
}
lines <- c("\\begin{tabular}{@{}rrr r rr rr@{}}", "\\toprule",
           tr("", "", "", "", "\\multicolumn{2}{c}{Permutation}", "\\multicolumn{2}{c}{Tuned $q$}"),
           "\\cmidrule(lr){5-6}\\cmidrule(l){7-8}",
           tr("$n$", "Censored", "$p$", "Events", "Empty", "All failed", "Empty", "All failed"), "\\midrule")
for (sc in SC) {
  k <- conf(sc); z <- sg[sg$scenario == sc, ]
  ref <- k$n == cfg$n_ref && k$cens == cfg$censor_ref && k$p == cfg$p_ref
  r <- sim_row(z, c(paste0(int(k$n), if (ref) "$^{\\dagger}$" else ""), sprintf("%s\\%%", pct(k$cens, 0)), int(k$p)))
  lines <- c(lines, tr(r))
}
emp <- vapply(c("survfcbf", "survfcbf_cvq"), function(s) mean(is.na(sg$sensitivity[sg$selector == s])), 0)
agree("empty-panel share", emp, unlist(js$simulation$empty_panel_share[c("survfcbf", "survfcbf_cvq")]))
agree("simulated datasets", sum(sg$selector == "survfcbf"), js$simulation$n_replicates, 0)
all_rows <- sim_row(sg, c("All", "", ""))
all_rows[4] <- ""
lines <- c(lines, "\\midrule", tr(all_rows), foot)
put("tab_simulation", lines)
mac("SuppSimConfigs", int(length(SC))); mac("SuppSimReps", int(cfg$R))
mac("SuppSimDatasets", int(sum(sg$selector == "survfcbf")))
mac("SuppSimEmptyPerm", pct(emp[["survfcbf"]]))

## The grid ran with ties broken by column order. grid_tiebreak_check.R reran the
## permutation mode with both rules at the same delta on part of it (Section S10).
tb <- utils::read.csv(file.path(RES, "grid_tiebreak_check.csv"))
stopifnot("the tie-break changed a grid selection" = all(tb$same_set) && all(tb$same_order),
          "a relevant feature tied in the grid check" = sum(tb$n_tied_relevant) == 0L)
tb_ds <- unique(tb[c("scenario", "replicate")])
mac("SuppTieCheckDatasets", int(nrow(tb_ds)))
mac("SuppTieCheckConfigs", int(length(unique(tb$scenario))))
mac("SuppTieCheckReps", int(length(unique(tb$replicate))))
mac("SuppTieCheckN", int(nrow(tb)))
mac("SuppTieCheckMulti", int(sum(tb$n_relevant >= 2L)))

## ---- Table: validation against the FCBF package ---------------------------------------
xv <- utils::read.csv(file.path(RES, "validation/fcbf_crossvalidation.csv"), check.names = FALSE)
sc_genes <- as.integer(sub(" genes.*", "", xv$detail[xv$check == "selections_scdengue"]))
sc_cells <- as.integer(sub(".* x ([0-9]+) cells.*", "\\1", xv$detail[xv$check == "selections_scdengue"]))
VAL <- list(
  su_random_pairs = "SU of random factor pairs ($n$ 20 to 120, 2 to 6 levels)",
  selections_nondegenerate = "Selections, 400 generated datasets, those without an exact SU tie",
  selections_with_ties = "Selections, 300 generated datasets with duplicated features allowed",
  selections_scdengue = sprintf("Selections on scDengue (%s genes, %d cells), five thresholds",
                                int(sc_genes), sc_cells),
  topk_random = "400 generated datasets, $K$ = 1 to 6, one of four thresholds each",
  topk_relevance_weights = "With precomputed relevance and observation weights",
  topk_ties = "Duplicated columns and tied relevance",
  topk_scdengue = "scDengue, threshold 0.05, $K$ = 1 to 10")
stopifnot("validation checks changed" = setequal(names(VAL), xv$check))
vrow <- function(ck) { z <- xv[xv$check == ck, ]
  tr(VAL[[ck]], int(z$comparisons), int(z$mismatches), int(z$mismatches_at_tie),
     if (ck == "selections_nondegenerate") int(z$skipped_ties) else "") }
topk <- grep("^topk_", xv$check, value = TRUE)
agree("top-K comparisons", sum(xv$comparisons[xv$check %in% topk]), js$validation$topk_comparisons, 0)
lines <- c("\\begin{tabularx}{\\linewidth}{@{}X rrrr@{}}", "\\toprule",
           tr("Check", "Compared", "Mismatches", "At a tie", "Skipped"), "\\midrule",
           "\\multicolumn{5}{@{}l}{\\textit{A. \\pkg{survfcbf} against the \\pkg{FCBF} package}} \\\\",
           vapply(c("su_random_pairs", "selections_nondegenerate", "selections_with_ties",
                    "selections_scdengue"), vrow, ""),
           "\\midrule",
           "\\multicolumn{5}{@{}l}{\\textit{B. top-$K$ run against the first $K$ features of the full run}} \\\\",
           vapply(c("topk_random", "topk_relevance_weights", "topk_ties", "topk_scdengue"), vrow, ""),
           "\\midrule", tr("All top-$K$ checks", int(sum(xv$comparisons[xv$check %in% topk])),
                           int(sum(xv$mismatches[xv$check %in% topk])), "", ""),
           "\\bottomrule", "\\end{tabularx}")
put("tab_validation", lines)
mac("SuppSuMaxDiff", sprintf("%s \\times 10^{%d}", num(signif(xv$max_abs_diff[xv$check == "su_random_pairs"], 3) /
      10^floor(log10(xv$max_abs_diff[xv$check == "su_random_pairs"])), 2),
      as.integer(floor(log10(xv$max_abs_diff[xv$check == "su_random_pairs"])))))
mac("SuppScGenes", int(sc_genes)); mac("SuppScCells", int(sc_cells))

## ---- Table: runtime -------------------------------------------------------------------
rt <- utils::read.csv(file.path(RES, "realdata_runtime/runtime_scaling.csv"), check.names = FALSE)
fit <- stats::lm(log(survfcbf_topk) ~ log(p), data = rt)
agree("runtime exponent", unname(coef(fit)[2]), js$runtime$scaling_exponent)
RT <- c(survfcbf_topk = "\\pkg{survfcbf} top-$K$", variance = "Variance", carss = "CAR", mim_surv = "MIM")
## real data: mean selection time per outer fold of top-K in the TCGA benchmark
rtk <- rd[rd$selector == "survfcbf_topk", ]
mac("SuppRealTimeMin", num(min(rtk$runtime), 1)); mac("SuppRealTimeMax", num(max(rtk$runtime), 1))
mac("SuppRealPMin", int(min(rtk$p))); mac("SuppRealPMax", int(max(rtk$p)))
mac("SuppMetabricTime", num(mb$runtime[mb$selector == "survfcbf_topk"], 1))
lines <- c("\\begin{tabular}{@{}r rrrr@{}}", "\\toprule", tr("$p$", RT), "\\midrule")
for (i in seq_len(nrow(rt))) lines <- c(lines, tr(int(rt$p[i]), num(unlist(rt[i, names(RT)]), 2)))
put("tab_runtime", c(lines, foot))
mac("SuppRuntimeExponent", num(coef(fit)[2], 2))

## ---- Table: software versions, from the lockfile ----------------------------------------
lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
SW <- list(
  c("survival", "Null Cox model, martingale residuals"),
  c("mlr3", "Learner and measure framework"),
  c("mlr3proba", "Survival tasks, Uno's C, integrated Brier score"),
  c("mlr3extralearners", "Cross-validated penalised Cox learner"),
  c("glmnet", "Ridge and lasso Cox models"),
  c("carSurv", "CAR scores"),
  c("praznik", "MIM and mRMR filters"),
  c("stabm", "Selection stability"),
  c("fgsea", "Gene-set enrichment"),
  c("msigdbr", "MSigDB gene sets"),
  c("org.Hs.eg.db", "Ensembl to symbol mapping"),
  c("FCBF", "Reference implementation for the validation"))
src_of <- function(p) {
  if (identical(p$Source, "GitHub"))
    sprintf("GitHub %s/%s@%s", p$RemoteUsername, p$RemoteRepo, substr(p$RemoteSha, 1, 7))
  else if (identical(p$Source, "Bioconductor")) "Bioconductor"
  else if (!is.null(p$Repository)) p$Repository else p$Source
}
lines <- c("\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}l l l l@{}}", "\\toprule",
           tr("Package", "Version", "Source", "Role"), "\\midrule")
for (s in SW) {
  p <- lock$Packages[[s[1]]]
  stopifnot("package missing from renv.lock" = !is.null(p))
  lines <- c(lines, tr(sprintf("\\pkg{%s}", s[1]), p$Version, src_of(p), s[2]))
}
put("tab_software", c(lines, foot))
mac("SuppLockR", lock$R$Version); mac("SuppLockBioc", lock$Bioconductor$Version)

## ---- Tables: redundancy, mRMR and q at matched size (realdata/panel_redundancy.R) --------
RR <- file.path(RES, "realdata_redundancy")
rr <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+\\.csv$", full.names = TRUE),
                            utils::read.csv, check.names = FALSE))
rr$K <- as.integer(sub(".*_K", "", rr$selector)); rr$family <- sub("_K[0-9]+$", "", rr$selector)
stopifnot("redundancy runs must cover every cohort" = setequal(unique(rr$cohort), COH),
          "redundancy runs lost folds" = all(rr$n_failed == 0))
RED <- c(survfcbf_q2 = "\\pkg{survfcbf}", mrmr = "mRMR", mrmr_q2 = "mRMR, $q = 2$", mim_surv_q2 = "MIM, $q = 2$",
         variance = "Variance", carss = "CAR", mim_surv = "MIM")
## the untruncated C of the shared rows reproduces the cardinality curve (checked in
## paper/make_statistics.R); the table gives the truncated C of the same single split
red_row <- function(z, lab, K) tr(lab, if (is.null(K)) "" else K, num(mean(z$uno_c_trunc)), num(mean(z$stability)),
                                  num(mean(z$mean_abs_r)), pct(mean(z$share_r05)), int(sum(z$n_null_fit)))
lines <- c("\\begin{tabular}{@{}l r rrrrr@{}}", "\\toprule",
           tr("Selector", "$K$", "Uno's C", "Stability", "Mean $|r|$", "$|r| \\ge 0.5$ (\\%)", "Null fits"),
           "\\midrule")
for (K in c(5L, 10L, 20L)) {
  for (f in names(RED)) lines <- c(lines, red_row(rr[rr$family == f & rr$K == K, ], RED[[f]], K))
  if (K == 10L) for (q in c(4L, 10L))
    lines <- c(lines, red_row(rr[rr$family == sprintf("survfcbf_q%d", q) & rr$K == K, ],
                              sprintf("\\pkg{survfcbf}, $q = %d$", q), K))
  if (K != 20L) lines <- c(lines, "\\addlinespace")
}
put("tab_redundancy", c(lines, foot))
k10 <- rr[rr$K == 10L, ]
mr <- function(f, col = "mean_abs_r") mean(k10[[col]][k10$family == f])
for (f in c("survfcbf_q2", "mim_surv", "mrmr", "carss", "mim_surv_q2", "mrmr_q2"))
  agree(sprintf("in-panel |r| at K=10, %s", f), mr(f), js$redundancy$mean_abs_r_K10[[f]])
agree("mRMR Uno's C at K=10", mr("mrmr", "uno_c"), js$redundancy$uno_c_K10$mrmr)
at10 <- function(f, col) { z <- k10[k10$family == f, ]; z[[col]][match(COH, z$cohort)] }
below <- sum(at10("survfcbf_q2", "mean_abs_r") < at10("mim_surv", "mean_abs_r"))
below2 <- sum(at10("survfcbf_q2", "mean_abs_r") < at10("mim_surv_q2", "mean_abs_r"))
belowm2 <- sum(at10("survfcbf_q2", "mean_abs_r") < at10("mrmr_q2", "mean_abs_r"))
agree("cohorts below MIM", below, js$redundancy$cohorts_survfcbf_below_mim, 0)
agree("cohorts below MIM q=2", below2, js$redundancy$cohorts_survfcbf_below_mim_q2, 0)
mac("SuppRedTopk", num(mr("survfcbf_q2"))); mac("SuppRedMim", num(mr("mim_surv")))
mac("SuppRedMrmr", num(mr("mrmr"))); mac("SuppRedCar", num(mr("carss")))
mac("SuppRedMimQtwo", num(mr("mim_surv_q2"))); mac("SuppRedMrmrQtwo", num(mr("mrmr_q2")))
mac("SuppRedBelowMim", word_of(below)); mac("SuppRedBelowMimQtwo", word_of(below2))
mac("SuppRedBelowMrmrQtwo", word_of(belowm2))
mac("SuppShareTopk", pct(mr("survfcbf_q2", "share_r05"))); mac("SuppShareMim", pct(mr("mim_surv_q2", "share_r05")))
## What the Markov-blanket trace of the ten-gene panels records (realdata/panel_groups.R)
PG <- file.path(RES, "realdata_panel_groups")
pg_rm <- do.call(rbind, lapply(sort(list.files(PG, "_removed\\.csv$", full.names = TRUE)), utils::read.csv))
pg_pn <- do.call(rbind, lapply(sort(list.files(PG, "_panel\\.csv$", full.names = TRUE)), utils::read.csv))
stopifnot("panel_groups.R did not cover the eleven cohorts" = length(unique(pg_pn$cohort)) == 11L,
          "a trace panel differs from the published panel" = all(pg_pn$panel_matches))
pg_co <- sort(unique(pg_pn$cohort))
pg_by <- function(fun, d) vapply(pg_co, function(co) fun(d[d$cohort == co, ]), 0)
pg_medr <- pg_by(function(z) stats::median(z$abs_r), pg_rm)
pg_bg   <- pg_by(function(z) mean(z$bg_abs_r), pg_pn)
pg_r05  <- pg_by(function(z) mean(z$abs_r >= 0.5), pg_rm)
pg_pctl <- pg_by(function(z) stats::median(z$pct_in_blanket), pg_rm)
pg_cls  <- pg_by(function(z) mean(z$blanket_is_closest), pg_rm)
stopifnot("pruned genes not above the background" = all(pg_medr > pg_bg))
mac("SuppGrpFolds", int(length(unique(paste(pg_pn$cohort, pg_pn$fold)))))
mac("SuppGrpRMin", num(min(pg_medr), 2)); mac("SuppGrpRMax", num(max(pg_medr), 2))
mac("SuppGrpBgMin", num(min(pg_bg), 2)); mac("SuppGrpBgMax", num(max(pg_bg), 2))
mac("SuppGrpHalfMin", pct(min(pg_r05), 0)); mac("SuppGrpHalfMax", pct(max(pg_r05), 0))
mac("SuppGrpTopPct", int(ceiling(100 * (1 - min(pg_pctl)))))
mac("SuppGrpClosestMin", pct(min(pg_cls), 0)); mac("SuppGrpClosestMax", pct(max(pg_cls), 0))
mac("SuppGrpEmpty", pct(mean(pg_pn$group_size == 0), 0))
mac("SuppGrpMaxGroup", int(max(pg_pn$group_size)))
## q against the default q = 2 at K = 10, per cohort, truncated C of the same split
qd <- function(q, col) at10(sprintf("survfcbf_q%d", q), col) - at10("survfcbf_q2", col)
mac("SuppQfourC", num(mean(qd(4L, "uno_c_trunc")))); mac("SuppQtenC", num(mean(qd(10L, "uno_c_trunc"))))
mac("SuppQdiffMin", num(min(qd(4L, "uno_c_trunc"), qd(10L, "uno_c_trunc"))))
mac("SuppQdiffMax", num(max(qd(4L, "uno_c_trunc"), qd(10L, "uno_c_trunc"))))
mac("SuppQfourLessStable", word_of(sum(qd(4L, "stability") < 0)))
mac("SuppQtenLessStable", word_of(sum(qd(10L, "stability") < 0)))
mac("SuppQtenMoreRedundant", word_of(sum(qd(10L, "mean_abs_r") > 0)))

## ---- Table: Cox on censoring (S2), against 100 random panels per cohort ----------------------
rv <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+_reverse\\.csv$", full.names = TRUE),
                            utils::read.csv, check.names = FALSE))
rnd <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+_reverse_random\\.csv$", full.names = TRUE),
                             utils::read.csv, check.names = FALSE))
stopifnot("censoring runs lost folds" = all(rv$n_failed == 0) && all(rnd$n_failed == 0),
          "100 random panels per cohort" = all(table(rnd$cohort) == 100L))
rv <- rv[rv$selector != "random_K10", ]
rv$pct <- mapply(function(co, v) mean(rnd$harrell_c_censoring[rnd$cohort == co] < v), rv$cohort, rv$harrell_c_censoring)
mcc <- utils::read.csv(file.path(RES, "realdata_metabric/censoring_cox.csv"), check.names = FALSE)
CENS <- c(survfcbf_q2_K10 = "\\pkg{survfcbf}, $q = 2$", survfcbf_q4_K10 = "\\pkg{survfcbf}, $q = 4$",
          survfcbf_q10_K10 = "\\pkg{survfcbf}, $q = 10$", mrmr_K10 = "mRMR", mrmr_q2_K10 = "mRMR, $q = 2$",
          mim_surv_K10 = "MIM", mim_surv_q2_K10 = "MIM, $q = 2$", variance_K10 = "Variance", carss_K10 = "CAR")
lines <- c("\\begin{tabular}{@{}l rrr rr r r@{}}", "\\toprule",
           tr("", "\\multicolumn{5}{c}{TCGA, eleven cohorts}", "", ""), "\\cmidrule(lr){2-6}",
           tr("Panel of ten genes", "Mean", "Min", "Max", "Percentile", "Above 95th", "METABRIC", "Null fits"), "\\midrule")
for (s in names(CENS)) {
  z <- rv[rv$selector == s, ]
  mp_ <- mcc$percentile_among_random[mcc$selector == s]
  lines <- c(lines, tr(CENS[[s]], num(mean(z$harrell_c_censoring)), num(min(z$harrell_c_censoring)),
                       num(max(z$harrell_c_censoring)), num(mean(z$pct), 2), int(sum(z$pct > 0.95)),
                       if (length(mp_)) num(mp_, 2) else "n.r.", int(sum(z$n_null_fit))))
}
rmc <- tapply(rnd$harrell_c_censoring, rnd$cohort, mean)
lines <- c(lines, "\\midrule", tr("100 random panels", num(mean(rnd$harrell_c_censoring)), num(min(rmc)), num(max(rmc)),
                                  "", "", "", int(sum(rnd$n_null_fit))))
put("tab_censoring", c(lines, foot))
cm <- tapply(rv$harrell_c_censoring, rv$selector, mean)
above <- tapply(rv$pct > 0.95, rv$selector, sum)
agree("censoring C, survfcbf", cm[["survfcbf_q2_K10"]], js$censoring_cox$harrell_c$survfcbf_q2_K10)
agree("random panels mean", mean(rnd$harrell_c_censoring), js$censoring_cox$random_mean)
agree("cohorts above p95", above[names(CENS)], unlist(js$censoring_cox$cohorts_above_random_p95[names(CENS)]), 0)
q2p <- rv[rv$selector == "survfcbf_q2_K10", ]
mac("SuppCensTopk", num(cm[["survfcbf_q2_K10"]])); mac("SuppCensRandomMean", num(mean(rnd$harrell_c_censoring)))
mac("SuppCensTopkMaxPct", num(max(q2p$pct), 2)); mac("SuppCensTopkMaxCohort", q2p$cohort[which.max(q2p$pct)])
mac("SuppCensAboveQfour", word_of(above[["survfcbf_q4_K10"]])); mac("SuppCensAboveQten", word_of(above[["survfcbf_q10_K10"]]))
mac("SuppCensAboveMim", word_of(above[["mim_surv_K10"]])); mac("SuppCensAboveMrmr", word_of(above[["mrmr_K10"]]))
ord <- function(p) { k <- round(100 * p); paste0(k, if (k %% 100 %in% 11:13) "th" else c("th", "st", "nd", "rd", rep("th", 6))[k %% 10 + 1]) }
mac("SuppCensMetabricTopk", ord(mcc$percentile_among_random[mcc$selector == "survfcbf_q2_K10"]))
mac("SuppCensMetabricQten", ord(mcc$percentile_among_random[mcc$selector == "survfcbf_q10_K10"]))
mac("SuppCensMetabricMrmr", ord(mcc$percentile_among_random[mcc$selector == "mrmr_K10"]))
stopifnot("the q = 2 panel exceeds the random 95th percentile; the note says it never does" =
            above[["survfcbf_q2_K10"]] == 0 && mcc$percentile_among_random[mcc$selector == "survfcbf_q2_K10"] <= 0.95)

## ---- Survival target, centre and follow-up (S1, S2) -----------------------------------------
ts <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_target_sensitivity"), pattern = "^[A-Z]+\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
tm <- function(col, tgt) mean(ts[[col]][ts$target == tgt])
agree("IPCW truncated C", c(tm("uno_c_trunc", "ipcw"), tm("uno_c_trunc", "martingale")),
      c(js$target$means$uno_c_trunc$ipcw, js$target$means$uno_c_trunc$martingale))
mac("SuppIpcwC", num(tm("uno_c_trunc", "ipcw"))); mac("SuppMartC", num(tm("uno_c_trunc", "martingale")))
mac("SuppIpcwCens", num(tm("harrell_c_censoring", "ipcw"))); mac("SuppMartCens", num(tm("harrell_c_censoring", "martingale")))
mac("SuppIpcwJaccardMax", num(max(ts$jaccard_between_targets), 0))
mac("SuppGradMart", num(js$target$gradient_spearman$martingale, 2)); mac("SuppGradIpcw", num(js$target$gradient_spearman$ipcw, 2))
mac("SuppBrcaMartCens", num(js$target$brca_censoring_c$martingale)); mac("SuppBrcaIpcwCens", num(js$target$brca_censoring_c$ipcw))
ce <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_center"), pattern = "^[A-Z]+\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
agree("centre C on censoring", mean(ce$c_centre_censoring), js$centre$c_centre_censoring)
mac("SuppCentreMin", int(min(ce$n_centres))); mac("SuppCentreMax", int(max(ce$n_centres)))
mac("SuppCentreCens", num(mean(ce$c_centre_censoring)))
mac("SuppCentreCensMin", num(min(ce$c_centre_censoring))); mac("SuppCentreCensMax", num(max(ce$c_centre_censoring)))
mac("SuppCentreDeath", num(mean(ce$c_centre_death)))
mac("SuppCentreDiffer", word_of(sum(ce$kruskal_p_followup_censored < 0.001)))
mac("SuppFollowMin", int(min(ce$median_followup_reverse_km))); mac("SuppFollowMax", int(max(ce$median_followup_reverse_km)))
fo <- do.call(rbind, lapply(list.files(RR, pattern = "^[A-Z]+_folds\\.csv$", full.names = TRUE),
                            utils::read.csv, check.names = FALSE))
agree("untruncated off by 0.10", mean(abs(fo$uno_c - fo$uno_c_trunc) > 0.10), js$repeated_cv$untruncated_off$share_gap_over_0.10)
mac("SuppTruncFolds", int(nrow(fo)))
mac("SuppTruncOffTen", pct(mean(abs(fo$uno_c - fo$uno_c_trunc) > 0.10)))
mac("SuppTruncOffTwenty", pct(mean(abs(fo$uno_c - fo$uno_c_trunc) > 0.20)))

## ---- Table: stability threshold (S5) ----------------------------------------------------------
th <- do.call(rbind, lapply(list.files(file.path(RES, "realdata_stability_theta"), pattern = "\\.csv$",
                                       full.names = TRUE), utils::read.csv, check.names = FALSE))
th$family <- sub("_K10$", "", th$selector)
THS <- sort(unique(th$theta), decreasing = TRUE)
lines <- c("\\begin{tabular}{@{}l rrr@{}}", "\\toprule",
           tr("Selector, $K = 10$", sprintf("$\\theta = %s$", THS)), "\\midrule")
for (f in names(RED)) lines <- c(lines, tr(RED[[f]], vapply(THS, function(t) num(mean(th$stability[th$family == f & th$theta == t])), "")))
gap_th <- function(t) { a <- th[th$family == "survfcbf_q2" & th$theta == t, ]; b <- th[th$family == "mim_surv_q2" & th$theta == t, ]
  m <- merge(a, b, by = "cohort"); c(mean(m$stability.x - m$stability.y), sum(m$stability.x < m$stability.y)) }
lines <- c(lines, "\\midrule", tr("\\pkg{survfcbf} minus MIM, $q = 2$", vapply(THS, function(t) num(gap_th(t)[1]), "")),
           tr("Cohorts with \\pkg{survfcbf} lower", vapply(THS, function(t) int(gap_th(t)[2]), "")), foot)
put("tab_theta", lines)
for (t in THS) agree(sprintf("theta %s gap", t), gap_th(t)[1], js$stability_theta$vs_mim_surv_q2[[as.character(t)]]$diff)
mac("SuppThetaGapNine", num(gap_th(0.9)[1])); mac("SuppThetaGapSeven", num(gap_th(0.7)[1]))
mac("SuppThetaGapFive", num(gap_th(0.5)[1]))
mac("SuppThetaLowerMin", word_of(min(gap_th(0.9)[2], gap_th(0.7)[2], gap_th(0.5)[2])))
mac("SuppThetaLowerMax", word_of(max(gap_th(0.9)[2], gap_th(0.7)[2], gap_th(0.5)[2])))
mac("SuppThetaDrawsLow", int(unique(th$draws[th$theta == min(THS)])))
mac("SuppThetaDraws", int(unique(th$draws[th$theta == max(THS)])))

## ---- Table: simulation at matched size (S7) ---------------------------------------------
sm <- utils::read.csv(file.path(RES, "raw_sim_matched_size.csv"), check.names = FALSE)
sm$K <- as.integer(sub(".*_K", "", sm$selector)); sm$family <- sub("_K[0-9]+$", "", sm$selector)
stopifnot("matched-size simulation must hold the grid's datasets" =
            nrow(unique(sm[, c("scenario", "replicate")])) == sum(sg$selector == "survfcbf"),
          "matched-size simulation lost folds" = all(sm$n_failed == 0))
mcse <- function(x) sprintf("%s (%s)", num(mean(x)), num(stats::sd(x) / sqrt(length(x))))
SIMCOL <- c("sensitivity", "fdr", "redundant_copies", "exact", "uno_c_trunc", "stability")
mac("SuppSimMrmrQtwoThree", num(mean(sm$sensitivity[sm$family == "mrmr_q2" & sm$K == 3L]), 2))
mac("SuppSimMrmrQtwoTen", num(mean(sm$sensitivity[sm$family == "mrmr_q2" & sm$K == 10L]), 2))
mac("SuppSimMrmrQtwoCopies", num(mean(sm$redundant_copies[sm$family == "mrmr_q2" & sm$K == 10L]), 2))
mac("SuppSimTopkThree", num(mean(sm$sensitivity[sm$family == "survfcbf_q2" & sm$K == 3L]), 2))
mac("SuppSimTopkTen", num(mean(sm$sensitivity[sm$family == "survfcbf_q2" & sm$K == 10L]), 2))
lines <- c("\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}l r rrrrrr@{}}", "\\toprule",
           tr("Selector", "$K$", "Recovered", "FDR", "Copies", "Exact", "Uno's C", "Stability"), "\\midrule")
for (K in c(3L, 10L)) {
  for (f in names(RED)) {
    z <- sm[sm$family == f & sm$K == K, ]
    lines <- c(lines, tr(RED[[f]], K, vapply(SIMCOL, function(v) mcse(z[[v]]), "")))
  }
  if (K == 3L) lines <- c(lines, "\\addlinespace")
}
put("tab_sim_matched", c(lines, foot))
for (K in c("3", "10"))
  agree(sprintf("sim recovery K=%s, survfcbf", K),
        mean(sm$sensitivity[sm$family == "survfcbf_q2" & sm$K == as.integer(K)]),
        js$sim_matched$sensitivity[[K]]$survfcbf_q2)
nf <- tapply(sm$n_null_fit, sm$family, sum)
mac("SuppSimMNullVar", int(nf[["variance"]])); mac("SuppSimMNullTopk", int(nf[["survfcbf_q2"]]))
mac("SuppSimMFolds", int(5L * nrow(unique(sm[, c("scenario", "replicate")])) * 2L))

## ---- claims made in words in supplement_body.tex ----------------------------------------
rsK <- do.call(rbind, lapply(names(MKD), function(K) { z <- utils::read.csv(file.path(RES, MKD[[K]], "rope_sensitivity.csv"))
  z$rival <- sub("_K[0-9]+$", "", z$rival); cbind(K = as.integer(K), z) }))
vh <- function(r, rope, m = "uno_c_trunc") rsh$verdict[rsh$rival == r & rsh$metric == m & rsh$rope == rope & rsh$summary == "argmax"]
stopifnot(
  "finer binnings must exceed the random 95th percentile in the same number of cohorts" =
    above[["survfcbf_q4_K10"]] == above[["survfcbf_q10_K10"]],
  "headline decided at ROPE 0.01 and 0.02, not at 0.03" =
    all(vapply(RK, function(r) vh(r, 0.01) == "rival better" && vh(r, 0.02) == "rival better" &&
                 vh(r, 0.03) != "rival better", TRUE)),
  "headline undecided with the untruncated C" =
    all(hinf$verdict[hinf$metric == "uno_c"] == "undecided at n=11"),
  "mRMR equivalent at ROPE 0.02 and 0.03 at every matched size" =
    all(rsK$verdict[rsK$rival %in% c("mrmr", "mrmr_q2") & rsK$metric == "uno_c_trunc" & rsK$rope %in% c(0.02, 0.03) &
                      rsK$summary == "argmax"] == "practically equivalent"),
  "the only decided matched C comparison is equivalence with CAR at K = 20" =
    identical(minf$K[minf$metric == "uno_c_trunc" & minf$rival %in% c("variance", "carss", "mim_surv") &
                       minf$verdict != "undecided at n=11"], 20L),
  "METABRIC gap in all ten splits" = all(vapply(mg, function(g) g$splits_topk_worse == 10L, TRUE)))

## ---- macros -----------------------------------------------------------------------------
writeLines(c("% generated by paper/make_supplementary_tables.R, do not edit",
             sprintf("\\newcommand{\\%s}{%s}", names(macros), unlist(macros))),
           file.path(OUT, "supp_numbers.tex"))
cat(sprintf("wrote %d tables and %d macros to %s\n",
            length(list.files(OUT, pattern = "^tab_")), length(macros), OUT))
cat("SUPPLEMENTARY_TABLES_DONE\n")
