# What each gene of a top-K panel stands in for, on the eleven TCGA cohorts. With the
# panel, survfcbf returns every gene it pruned and the panel gene that removed it as an
# approximate Markov blanket (sel$trace$removed). At q = 2 the relevance SU(Fj, Y) of a
# gene is small, so the blanket condition SU(Fi, Fj) >= SU(Fj, Y) can hold at modest
# correlation, and this measures how tight the resulting groups are. Same coding-only
# feature space, five outer folds and seed as panel_redundancy.R, and the panel must
# reproduce its survfcbf_q2_K10 panel in every fold (column panel_matches). Per fold, on
# the training rows:
#   <cohort>_removed.csv, one row per pruned gene: its blanket, the two SU values, the
#     absolute Pearson correlation with the blanket (abs_r), the percentile of that |r|
#     among the blanket's |r| with every other coding gene (pct_in_blanket), and whether
#     the blanket is the panel gene most correlated with it (blanket_is_closest);
#   <cohort>_panel.csv, one row per panel gene: its group size (genes it removed before
#     the cap of ten was reached) and its mean |r| with every other coding gene (bg_abs_r).
# One task per cohort. The published run (CESGA, 2026-11-20, 10 to 15 min per cohort, most
# of it reading the ARFF files, under 2 GB), with the reference panels of the v3 run:
#   RD_REF_DIR=results/realdata_redundancy_v3 sbatch -p short --qos=short --time=01:00:00 \
#     --mem=16G --cpus-per-task=1 --array=1-11 cesga/submit_script.sbatch realdata/panel_groups.R
# It reproduced the published panel in all 55 folds. Read by paper/make_supplementary_tables.R
# (Section S5).
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(farff); library(jsonlite); library(curl)
})
stopifnot("this needs the tie-break build of survfcbf" = "tiebreak" %in% names(formals(survfcbf_select)))
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(getwd(), "realdata/Data"))
co <- names(BOMMERT_IDS)[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))]
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_panel_groups")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

d0 <- load_bommert(co); X <- d0$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
task <- surv_task(d0$time, d0$status, X)
d <- survfcbf:::extract_surv(task)
folds <- survfcbf:::stratified_folds(d$status, 5L, 1L)
ref <- utils::read.csv(file.path(Sys.getenv("RD_REF_DIR", "results/realdata_redundancy"), paste0(co, "_panels.csv")))
ref <- ref[ref$selector == "survfcbf_q2_K10", ]
cat(sprintf("%s: n=%d events=%d p=%d\n", co, length(d$time), sum(d$status), ncol(d$X)))

removed <- list(); panel <- list()
for (f in seq_len(5L)) {
  tr <- which(folds != f)
  Xt <- as.matrix(d$X[tr, , drop = FALSE])
  fs <- as.integer((1 * 1000 + f) %% 2147483647)           # the fold seed of nested_cv_select
  s <- survfcbf:::with_seed(fs, suppressWarnings(
    survfcbf_select(d$time[tr], d$status[tr], Xt, q = 2L, top = 10L)))
  sel <- s$selected; feats <- colnames(Xt)
  refp <- ref$feature[ref$fold == f][order(ref$position[ref$fold == f])]
  rm_ <- s$trace$removed
  stopifnot(all(rm_$blanket %in% sel))
  R <- abs(stats::cor(Xt[, sel, drop = FALSE], Xt))         # panel genes x all coding genes
  R[cbind(seq_along(sel), sel)] <- NA                         # drop self-correlation
  if (nrow(rm_)) {
    bi <- match(rm_$blanket, sel)
    ar <- R[cbind(bi, rm_$feature)]
    removed[[f]] <- data.frame(
      cohort = co, fold = f, feature = feats[rm_$feature], blanket = feats[rm_$blanket],
      su_blanket_feature = rm_$su_blanket_feature, su_feature_target = rm_$su_feature_target,
      abs_r = ar,
      pct_in_blanket = vapply(seq_along(ar), function(i) mean(R[bi[i], ] <= ar[i], na.rm = TRUE), 0),
      blanket_is_closest = vapply(seq_along(ar), function(i) which.max(R[, rm_$feature[i]]) == bi[i], NA))
  }
  panel[[f]] <- data.frame(
    cohort = co, fold = f, position = seq_along(sel), feature = feats[sel],
    group_size = vapply(sel, function(g) sum(rm_$blanket == g), 0L),
    bg_abs_r = rowMeans(R, na.rm = TRUE),
    panel_matches = identical(feats[sel], refp))
}
removed <- do.call(rbind, removed); panel <- do.call(rbind, panel)
utils::write.csv(removed, file.path(outdir, paste0(co, "_removed.csv")), row.names = FALSE)
utils::write.csv(panel, file.path(outdir, paste0(co, "_panel.csv")), row.names = FALSE)
cat(sprintf("%s: panels match %d/5, pruned %d, median |r| with blanket %.3f, background %.3f\n",
            co, sum(tapply(panel$panel_matches, panel$fold, all)), NROW(removed),
            if (NROW(removed)) stats::median(removed$abs_r) else NA, mean(panel$bg_abs_r)))
utils::capture.output(utils::sessionInfo(), file = file.path(outdir, paste0(co, "_sessionInfo.txt")))
cat("PANEL_GROUPS_DONE", co, "\n")
