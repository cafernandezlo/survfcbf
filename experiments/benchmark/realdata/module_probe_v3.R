# GO/NO-GO probe #3 (best honest module-FCBF): stack CAR relevance onto modules and save
# the selected modules for enrichment. Fixed complete-linkage modules; per subsample rank
# them by eigengene SU relevance AND by eigengene carss-CAR relevance (the ablation showed
# CAR ~2x more stable than SU at gene level -- does it stack onto the module gain?).
# Also saves the FULL-cohort top-K module gene sets (Ensembl) so ORA/GSEA can test the
# interpretability pillar where the gene-level top-10 failed (0/11). SLURM array, coding-only.
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl); library(stabm) })
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")); co <- cohorts[id]
outdir <- "results/realdata_modules_v3"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".rds")); if (file.exists(out)) { cat(co, "skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
ens <- sub("_rna$", "", colnames(X))
cat(sprintf("%s: n=%d ev=%d p=%d\n", co, length(d$time), sum(d$status), ncol(X)))

hc  <- stats::hclust(stats::as.dist(1 - abs(stats::cor(X))), method = "complete")
mod <- stats::cutree(hc, k = 400L)
big <- as.integer(names(which(table(mod) >= 5L)))
mod_idx <- lapply(big, function(m) which(mod == m))
K <- 10L; B <- 20L
su_target  <- function(t, s) discretize_quantile(survival_target(t, s, method = "martingale"), q = 2L)$binned
eigengenes <- function(Xs) vapply(mod_idx, function(g) prcomp(Xs[, g, drop = FALSE], scale. = FALSE)$x[, 1], numeric(nrow(Xs)))

set.seed(1)
Lsu <- Lcar <- Lgene <- Lvar <- vector("list", B)
for (b in seq_len(B)) {
  i <- sample(length(d$time), round(0.8 * length(d$time)))
  Xs <- X[i, , drop = FALSE]; ts <- d$time[i]; ss <- d$status[i]
  EG <- eigengenes(Xs); yb <- su_target(ts, ss)
  rel_su  <- apply(EG, 2L, function(e) symmetrical_uncertainty(discretize_quantile(e, q = 2L)$binned, yb))
  rel_car <- abs(run_baseline(ent("carss"), ts, ss, EG, prop = NULL, seed = 1L)$score)
  Lsu[[b]]   <- order(rel_su,  decreasing = TRUE)[seq_len(K)]
  Lcar[[b]]  <- order(rel_car, decreasing = TRUE)[seq_len(K)]
  Lgene[[b]] <- survfcbf_select(ts, ss, Xs, q = 2L, top = K)$selected
  Lvar[[b]]  <- cut_ranking(apply(Xs, 2L, stats::var), k = K)
}
## full-cohort selected modules (for enrichment): gene Ensembl per selected module
EGf <- eigengenes(X); ybf <- su_target(d$time, d$status)
rel_su_f  <- apply(EGf, 2L, function(e) symmetrical_uncertainty(discretize_quantile(e, q = 2L)$binned, ybf))
rel_car_f <- abs(run_baseline(ent("carss"), d$time, d$status, EGf, prop = NULL, seed = 1L)$score)
top_su  <- order(rel_su_f,  decreasing = TRUE)[seq_len(K)]
top_car <- order(rel_car_f, decreasing = TRUE)[seq_len(K)]
mods_su  <- lapply(top_su,  function(m) ens[mod_idx[[m]]])
mods_car <- lapply(top_car, function(m) ens[mod_idx[[m]]])

res <- list(
  summary = data.frame(cohort = co, n = length(d$time), events = sum(d$status), p = ncol(X),
    n_modules = length(big),
    stab_module_su  = stabm::stabilityNogueira(Lsu,   p = length(big)),
    stab_module_car = stabm::stabilityNogueira(Lcar,  p = length(big)),
    stab_gene       = stabm::stabilityNogueira(Lgene, p = ncol(X)),
    stab_var        = stabm::stabilityNogueira(Lvar,  p = ncol(X)), stringsAsFactors = FALSE),
  mods_su = mods_su, mods_car = mods_car,
  background = ens)
saveRDS(res, out)
s <- res$summary
cat(sprintf("%s: module-SU %.3f | module-CAR %.3f | gene %.3f | variance %.3f\n",
            co, s$stab_module_su, s$stab_module_car, s$stab_gene, s$stab_var))
cat("MODULE_V3_DONE\n")
