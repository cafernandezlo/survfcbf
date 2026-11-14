# Follow-up and the tissue source site (TSS) in the TCGA cohorts. Ten random genes predict
# censoring above chance in every cohort (panel_redundancy.R), which suggests that follow-up
# varies with something that also shapes expression. The likeliest candidate is the
# centre that collected the sample, coded in the patient barcode (TCGA-<TSS>-<patient>).
# For each cohort this reports the follow-up (median by reverse Kaplan-Meier, range of the
# time variable, which is in days), whether follow-up of the censored patients differs
# between centres (Kruskal-Wallis), and the cross-validated Harrell's C of a Cox model on
# the centre alone, for death and for censoring, on the five stratified folds (seed 1) of
# the benchmark. Centres with fewer than five training patients are pooled, on the
# training set, into one level.
#
# SLURM array, one task per cohort, light (minutes):
#   sbatch -p short --qos=short --time=00:30:00 --mem=16G --cpus-per-task=1 --array=1-11 \
#     cesga/submit_script.sbatch realdata/center_followup.R
# Writes <cohort>.csv to RD_OUTDIR (default results/realdata_center).
suppressMessages({ library(survival); library(farff); library(jsonlite); library(curl) })
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(getwd(), "realdata/Data"))
co <- names(BOMMERT_IDS)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_center")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d <- load_bommert(co)
stopifnot("the loader must return the patient barcode" = length(d$barcode) == length(d$time))
tss <- substr(d$barcode, 6L, 7L)
time <- d$time; status <- d$status

## reverse Kaplan-Meier median follow-up
rkm <- survfit(Surv(time, 1L - status) ~ 1)
med_fu <- unname(summary(rkm)$table["median"])
cens <- status == 0L
kw <- if (length(unique(tss[cens])) > 1L) stats::kruskal.test(time[cens], factor(tss[cens]))$p.value else NA_real_

## cross-validated Harrell's C of the centre alone, on the benchmark's folds
folds <- survfcbf:::stratified_folds(status, 5L, 1L)
cv_c <- function(event) {
  vapply(1:5, function(f) {
    tr <- folds != f; te <- folds == f
    tab <- table(tss[tr]); keep <- names(tab)[tab >= 5L]
    lev <- function(x) factor(ifelse(x %in% keep, x, "pooled"), levels = c(keep, "pooled"))
    dtr <- data.frame(time = time[tr], ev = event[tr], g = droplevels(lev(tss[tr])))
    if (nlevels(dtr$g) < 2L) return(NA_real_)
    fit <- suppressWarnings(coxph(Surv(time, ev) ~ g, data = dtr))
    gte <- lev(tss[te]); gte[!gte %in% levels(dtr$g)] <- NA
    lp <- stats::predict(fit, newdata = data.frame(g = factor(gte, levels = levels(dtr$g))), type = "lp")
    lp[is.na(lp)] <- 0                      # a centre unseen in training gets the baseline risk
    survival::concordance(Surv(time[te], event[te]) ~ lp, reverse = TRUE)$concordance
  }, numeric(1))
}
c_death <- cv_c(status); c_cens <- cv_c(1L - status)

res <- data.frame(cohort = co, n = length(time), events = sum(status), n_centres = length(unique(tss)),
                  largest_centre_share = max(table(tss)) / length(tss),
                  time_min = min(time), time_median = stats::median(time), time_max = max(time),
                  median_followup_reverse_km = med_fu, kruskal_p_followup_censored = kw,
                  c_centre_death = mean(c_death), c_centre_death_se = stats::sd(c_death) / sqrt(5),
                  c_centre_censoring = mean(c_cens), c_centre_censoring_se = stats::sd(c_cens) / sqrt(5))
utils::write.csv(res, file.path(outdir, paste0(co, ".csv")), row.names = FALSE)
print(res, row.names = FALSE)
cat("CENTER_FOLLOWUP_DONE\n")
