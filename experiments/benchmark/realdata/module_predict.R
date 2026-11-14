# Module-FCBF prediction check: is the module-eigengene signature predictively
# within the CD of the baselines? Reuses the nested-CV harness by treating the module EIGENGENES as the feature matrix -- module
# selection is then just picking eigengene columns, which benchmark_selectors supports. Two
# variants: module-SU (SU of the eigengene vs the martingale target) and module-CAR (IPCW
# CAR of the eigengene). L2-Cox downstream, Uno's C. NOTE: for this quick check the modules
# (hence eigengenes) are fixed on the full cohort -- a mild leakage that makes the estimate
# slightly optimistic; module_leaksafe.R builds modules per training fold. SLURM array.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners)
  library(stabm); library(farff); library(jsonlite); library(curl)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")); co <- cohorts[id]
outdir <- "results/realdata_modules_v3"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, "_pred.csv")); if (file.exists(out)) { cat(co, "skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
hc  <- stats::hclust(stats::as.dist(1 - abs(stats::cor(X))), method = "complete")
mod <- stats::cutree(hc, k = 400L); big <- as.integer(names(which(table(mod) >= 5L)))
mod_idx <- lapply(big, function(m) which(mod == m))
EG <- vapply(mod_idx, function(g) prcomp(X[, g, drop = FALSE], scale. = FALSE)$x[, 1], numeric(nrow(X)))
colnames(EG) <- paste0("ME", seq_len(ncol(EG)))
cat(sprintf("%s: n=%d ev=%d p=%d -> %d module eigengenes\n", co, length(d$time), sum(d$status), ncol(X), ncol(EG)))

K <- 10L
su_target <- function(t, s) discretize_quantile(survival_target(t, s, method = "martingale"), q = 2L)$binned
sel_su  <- function(t, s, E) { yb <- su_target(t, s)
  order(apply(E, 2L, function(e) symmetrical_uncertainty(discretize_quantile(e, q = 2L)$binned, yb)), decreasing = TRUE)[seq_len(K)] }
sel_car <- function(t, s, E) order(abs(run_baseline(ent("carss"), t, s, E, prop = NULL, seed = 1L)$score), decreasing = TRUE)[seq_len(K)]

task <- surv_task(d$time, d$status, EG)
L <- mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")   # 1se -> null model (?score_benchmark)
tab <- benchmark_selectors(task, list(module_su = sel_su, module_car = sel_car),
                           L$clone(deep = TRUE), outer_folds = 5L, manifest = NULL, seed = 1L)
res <- data.frame(cohort = co, n = length(d$time), events = sum(d$status), p = ncol(X),
                  tab, row.names = NULL, check.names = FALSE)
utils::write.csv(res, out, row.names = FALSE)
cat(sprintf("%s: module_su uno_c %.3f | module_car uno_c %.3f\n",
            co, res$uno_c[res$selector == "module_su"], res$uno_c[res$selector == "module_car"]))
cat("MODULE_PRED_DONE\n")
