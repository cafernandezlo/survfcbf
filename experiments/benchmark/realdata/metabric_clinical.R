# Clinical anchor for the METABRIC external validation. Compares the Uno's C of a standard
# clinical prognostic model with the 10-gene survfcbf panel and with both combined, under
# the same leakage-safe harness as metabric_benchmark.R: nested_cv_select() with 5 outer
# folds, seed 1 and the default (sequential) rng, scored by score_benchmark(), ridge Cox
# learner surv.cv_glmnet (alpha = 0, s = "lambda.min").
#
# The three models are selector functions over ONE task that holds the clinical columns
# and the expression matrix, so they share the patients, the folds and the per-fold seeds,
# and their per-fold scores can be paired.
#   clinical             always returns the clinical columns.
#   panel                survfcbf_select(time, status, Xexpr, q = 2, top = 10)$selected,
#                        fitted on the expression columns of the training fold only.
#   clinical_plus_panel  the union of the two.
#
# Clinical variables, taken from metabric_clinical.rds (cBioPortal patient-level table).
# They are the standard breast-cancer prognostic set available in that table.
#   AGE_AT_DIAGNOSIS  age at diagnosis, years, numeric.
#   NPI               Nottingham prognostic index, numeric. It already combines tumour
#                     size, histological grade and nodal status, so those three are not
#                     added separately.
#   ER_IHC            ER status by immunohistochemistry, 1 = positive ("Positve" is the
#                     spelling as served), 0 = negative.
#   HER2_SNP6         HER2 status. The patient-level table has no IHC or FISH HER2 call,
#                     so this is the ERBB2 copy-number call from the SNP6 arrays,
#                     1 = GAIN, 0 = NEUTRAL or LOSS. UNDEF is treated as missing. GAIN is
#                     broader than clinical HER2 positivity (it includes low-level gains).
#   CLAUDIN_SUBTYPE   PAM50 subtype with claudin-low, one-hot with LumA (the largest and
#                     most favourable group) as reference: LumB, Her2, Basal, claudin-low,
#                     Normal. NC (not classified, a handful of patients) is treated as
#                     missing, not as a category. Normal is large enough to keep.
# PAM50 is itself expression-derived (a 50-gene classifier), so the clinical model already
# carries some transcriptomic information. Treatment variables (chemotherapy, hormone and
# radiotherapy) are left out on purpose: they are decisions taken after diagnosis and
# confounded by indication, not baseline prognostic factors.
#
# Patients: the METABRIC cohort of load_metabric() (OS endpoint, SHA-256 checked inputs)
# restricted to complete cases in the clinical variables. No imputation. The encoding is
# a fixed, outcome-free recoding with an a priori reference category, so nothing is learned
# from the data outside the training fold (glmnet standardises inside each training fold).
# Because of the complete-case restriction the folds differ from benchmark_metabric.csv,
# so the panel row here does not reproduce that file exactly.
#
# Output: results/realdata_metabric/clinical_reference.csv with one row per model plus the
# paired per-fold difference clinical_plus_panel minus clinical (mean over folds, with its
# standard error sd / sqrt(folds) in uno_c_se, and n_failed = folds that could not be
# paired). The console also prints a sensitivity with Uno's C and the IBS truncated at
# 120 months, because the untruncated IPCW measures degenerate in a fold whose test set
# holds follow-up beyond the training censoring support (see the comment at that block).
# METABRIC_DIR points at the .rds data, BENCH_DIR at experiments/benchmark.
suppressMessages({
  library(survfcbf)
  library(mlr3)
  library(mlr3proba)
  library(mlr3extralearners)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench)) setwd(bench)
source("realdata/load_metabric.R")

cfg <- list(outer = 5L, seed = 1L, q = 2L, topk = 10L)
OUTFILE <- "results/realdata_metabric/clinical_reference.csv"

t0 <- proc.time()[["elapsed"]]
d  <- load_metabric()
cat(sprintf("METABRIC (loader): n=%d events=%d p=%d\n", length(d$time), sum(d$status), ncol(d$X)))

## clinical table, same SHA-256 check as the loader
clin_file <- file.path(Sys.getenv("METABRIC_DIR"), "metabric_clinical.rds")
clin_sha  <- eval(formals(load_metabric)$clinical_sha256)
if (!identical(unname(tools::sha256sum(clin_file)), clin_sha))
  stop("METABRIC clinical SHA-256 mismatch: ", clin_file, call. = FALSE)
cl <- readRDS(clin_file)
cl <- cl[match(rownames(d$X), cl$patientId), ]
stopifnot("clinical rows must align with the expression rows" =
            identical(cl$patientId, rownames(d$X)))

## fixed recoding. Unexpected spellings stop the run instead of silently becoming NA.
stopifnot(
  "unexpected ER_IHC values"          = all(cl$ER_IHC %in% c("Positve", "Negative", NA)),
  "unexpected HER2_SNP6 values"       = all(cl$HER2_SNP6 %in% c("GAIN", "NEUTRAL", "LOSS", "UNDEF", NA)),
  "unexpected CLAUDIN_SUBTYPE values" = all(cl$CLAUDIN_SUBTYPE %in%
    c("LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal", "NC", NA))
)
recode01 <- function(x, one, zero) ifelse(x %in% one, 1, ifelse(x %in% zero, 0, NA_real_))
pam50_ref    <- "LumA"
pam50_levels <- c("LumB", "Her2", "Basal", "claudin-low", "Normal")
pam50 <- cl$CLAUDIN_SUBTYPE
cat("PAM50 counts in the loader cohort (NC is treated as missing):\n")
print(table(pam50, useNA = "ifany"))
pam50[pam50 == "NC"] <- NA
clin <- cbind(
  clin_age       = suppressWarnings(as.numeric(cl$AGE_AT_DIAGNOSIS)),
  clin_npi       = suppressWarnings(as.numeric(cl$NPI)),
  clin_er_pos    = recode01(cl$ER_IHC, "Positve", "Negative"),
  clin_her2_gain = recode01(cl$HER2_SNP6, "GAIN", c("NEUTRAL", "LOSS")),
  vapply(pam50_levels, function(l) ifelse(is.na(pam50), NA_real_, as.numeric(pam50 == l)),
         numeric(length(pam50)))
)
colnames(clin)[-(1:4)] <- paste0("clin_pam50_", gsub("-", "_", pam50_levels))
stopifnot("clinical column names collide with gene symbols" = !any(colnames(clin) %in% colnames(d$X)))

cat("Missing values per clinical column (loader cohort):\n")
print(colSums(is.na(clin)))
cc <- stats::complete.cases(clin)
time <- d$time[cc]
status <- d$status[cc]
n <- length(time)
events <- sum(status)
cat(sprintf("Complete cases: n=%d events=%d (dropped %d of %d)\n", n, events, sum(!cc), length(cc)))

X <- cbind(clin[cc, , drop = FALSE], d$X[cc, , drop = FALSE])
clin_names <- colnames(clin)
rm(d, cl)
invisible(gc())
task <- surv_task(time, status, X)
rm(X)
invisible(gc())

## selectors. mlr3 sorts the task's feature names, so the X a selector receives is not in
## the column order built above. Clinical and expression columns are therefore found by
## name, never by position. The panel sees only the expression columns of the training
## fold (survfcbf_select's Spearman tie-break makes it invariant to their order).
sel_clinical <- function(ti, st, X) match(clin_names, colnames(X))
sel_panel <- function(ti, st, X) {
  expr_idx <- which(!colnames(X) %in% clin_names)
  expr_idx[survfcbf_select(ti, st, X[, expr_idx, drop = FALSE], q = cfg$q, top = cfg$topk)$selected]
}
sel_both <- function(ti, st, X) union(sel_clinical(ti, st, X), sel_panel(ti, st, X))
models <- list(clinical = sel_clinical, panel = sel_panel, clinical_plus_panel = sel_both)

L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")   # 1se -> null model (?score_benchmark)
bms <- list()
scs <- list()
for (nm in names(models)) {
  t1 <- proc.time()[["elapsed"]]
  bms[[nm]] <- nested_cv_select(task, models[[nm]], L$clone(deep = TRUE),
                                outer_folds = cfg$outer, seed = cfg$seed)
  scs[[nm]] <- score_benchmark(bms[[nm]])
  cat(sprintf("%-20s done in %.1f min\n", nm, (proc.time()[["elapsed"]] - t1) / 60))
}

## each model must carry exactly what it claims: the clinical columns, a panel of genes
## only, and the union of the two from the same fold
for (f in seq_len(cfg$outer)) {
  fc <- bms$clinical$folds[[f]]$features
  fp <- bms$panel$folds[[f]]$features
  fb <- bms$clinical_plus_panel$folds[[f]]$features
  stopifnot("clinical model must use exactly the clinical columns" = setequal(fc, clin_names),
            "panel must contain expression columns only" = !any(fp %in% clin_names) && length(fp) == cfg$topk,
            "combined model must be clinical plus the panel of the same fold" = setequal(fb, union(fc, fp)))
}

## Uno's C truncated at the 80th percentile of each training fold's follow-up, the primary
## discrimination measure of the benchmark (realdata/uno_truncated.R)
source("realdata/uno_truncated.R")
ufs <- lapply(bms, uno_by_fold)
summ <- function(nm) {
  sc <- scs[[nm]]
  m  <- stats::setNames(sc$summary$mean, sc$summary$metric)
  se <- stats::setNames(sc$summary$se, sc$summary$metric)
  ut <- ufs[[nm]]$uno_c_trunc
  data.frame(model = nm, n = n, events = events, uno_c = m[["uno_c"]], uno_c_se = se[["uno_c"]],
             uno_c_trunc = mean(ut), uno_c_trunc_se = stats::sd(ut) / sqrt(length(ut)),
             harrell_c = m[["harrell_c"]], ibs = m[["ibs"]],
             n_failed = sc$n_failed, n_null_fit = sc$n_null_fit, row.names = NULL)
}
res <- do.call(rbind, lapply(names(models), summ))

## paired per-fold difference, clinical_plus_panel minus clinical (same folds, same seeds)
mets <- c("uno_c", "harrell_c", "ibs")
dd <- scs$clinical_plus_panel$per_fold[mets] - scs$clinical$per_fold[mets]
dt <- ufs$clinical_plus_panel$uno_c_trunc - ufs$clinical$uno_c_trunc
ok <- stats::complete.cases(dd)
diff_row <- data.frame(model = "clinical_plus_panel_minus_clinical", n = n, events = events,
                       uno_c = mean(dd$uno_c[ok]),
                       uno_c_se = stats::sd(dd$uno_c[ok]) / sqrt(sum(ok)),
                       uno_c_trunc = mean(dt), uno_c_trunc_se = stats::sd(dt) / sqrt(length(dt)),
                       harrell_c = mean(dd$harrell_c[ok]), ibs = mean(dd$ibs[ok]),
                       n_failed = sum(!ok), n_null_fit = NA_integer_, row.names = NULL)
res <- rbind(res, diff_row)

cat("\nPer-fold scores (Uno's C) and cardinality:\n")
per_fold <- data.frame(
  fold = seq_len(cfg$outer),
  k_clinical = scs$clinical$per_fold$cardinality,
  k_panel = scs$panel$per_fold$cardinality,
  k_both = scs$clinical_plus_panel$per_fold$cardinality,
  uno_clinical = scs$clinical$per_fold$uno_c,
  uno_panel = scs$panel$per_fold$uno_c,
  uno_both = scs$clinical_plus_panel$per_fold$uno_c,
  diff_both_minus_clinical = dd$uno_c,
  harrell_clinical = scs$clinical$per_fold$harrell_c,
  harrell_both = scs$clinical_plus_panel$per_fold$harrell_c
)
print(per_fold, digits = 4)
cat("\nPanel genes per fold:\n")
for (f in seq_len(cfg$outer))
  cat(sprintf("  fold %d: %s\n", f, paste(bms$panel$folds[[f]]$features, collapse = ", ")))
cat("\nResults:\n")
print(res, digits = 4)
cat(sprintf("\nUno's C, clinical_plus_panel minus clinical: %+.4f (paired SE %.4f, %d folds)\n",
            diff_row$uno_c, diff_row$uno_c_se, sum(ok)))

## Sensitivity, console only. The CSV keeps the harness measures of score_benchmark().
## The untruncated IPCW measures degenerate when a test fold holds follow-up beyond the
## last training censoring time: the training censoring KM reaches 0 there and a few pairs
## take almost all the weight (here fold 1 holds the longest follow-up of the cohort).
## Rescoring the same predictions truncated at 120 months, the usual 10-year horizon in
## breast cancer, fixed a priori, keeps the censoring weights away from 0.
TMAX <- 120
m_tmax <- list(uno_c = mlr3::msr("surv.cindex", weight_meth = "G2", t_max = TMAX, id = "uno_c"),
               ibs   = mlr3::msr("surv.graf", t_max = TMAX, id = "ibs"))
sc_t <- lapply(bms, score_benchmark, measures = m_tmax)
cat(sprintf("\nSensitivity, Uno's C and IBS truncated at t_max = %d months (not in the CSV):\n", TMAX))
print(data.frame(
  fold = seq_len(cfg$outer),
  uno_clinical = sc_t$clinical$per_fold$uno_c,
  uno_panel = sc_t$panel$per_fold$uno_c,
  uno_both = sc_t$clinical_plus_panel$per_fold$uno_c,
  diff_both_minus_clinical = sc_t$clinical_plus_panel$per_fold$uno_c - sc_t$clinical$per_fold$uno_c
), digits = 4)
print(do.call(rbind, lapply(names(sc_t), function(nm) {
  s <- sc_t[[nm]]$summary
  data.frame(model = nm, uno_c = s$mean[s$metric == "uno_c"], uno_c_se = s$se[s$metric == "uno_c"],
             ibs = s$mean[s$metric == "ibs"], ibs_se = s$se[s$metric == "ibs"], row.names = NULL)
})), digits = 4)
dt <- sc_t$clinical_plus_panel$per_fold$uno_c - sc_t$clinical$per_fold$uno_c
dt <- dt[is.finite(dt)]
cat(sprintf("Uno's C (t_max = %d), clinical_plus_panel minus clinical: %+.4f (paired SE %.4f, %d folds)\n",
            TMAX, mean(dt), stats::sd(dt) / sqrt(length(dt)), length(dt)))

dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, OUTFILE, row.names = FALSE)
cat(sprintf("METABRIC clinical reference done: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
cat("METABRIC_CLINICAL_DONE\n")
