# GO/NO-GO probe #2 (the acid test): does module-level selection stay stable at LOW events,
# where gene-level survfcbf collapsed to ~0.01-0.10? SLURM array over the 11 Bommert
# cohorts (coding-only, so the correlation matrix is ~16k not ~23k). Modules fixed on the
# full cohort (complete linkage, balanced, no blob), re-ranked per subsample by eigengene
# SU relevance. Compares module-FCBF vs gene survfcbf vs variance stability (Nogueira,
# top-10) per cohort. If module-FCBF holds ~>=0.5 across the low-event cohorts, the pivot
# is validated on the regime that broke the gene-level method.
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl); library(stabm) })
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")); co <- cohorts[id]
outdir <- "results/realdata_modules"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".rds")); if (file.exists(out)) { cat(co, "skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
cat(sprintf("%s: n=%d ev=%d p=%d (coding)\n", co, length(d$time), sum(d$status), ncol(X)))

## modules: complete linkage on 1-|cor| (balanced, no blob), fixed on full cohort
hc  <- stats::hclust(stats::as.dist(1 - abs(stats::cor(X))), method = "complete")
mod <- stats::cutree(hc, k = 400L)
big <- as.integer(names(which(table(mod) >= 5L)))
mod_idx <- lapply(big, function(m) which(mod == m))
K <- 10L; B <- 20L
su_target <- function(t, s) discretize_quantile(survival_target(t, s, method = "martingale"), q = 2L)$binned
sel_module <- function(t, s, Xs) {
  yb <- su_target(t, s)
  rel <- vapply(mod_idx, function(g) symmetrical_uncertainty(
    discretize_quantile(prcomp(Xs[, g, drop = FALSE], scale. = FALSE)$x[, 1], q = 2L)$binned, yb), numeric(1))
  order(rel, decreasing = TRUE)[seq_len(K)]
}
set.seed(1)
Lmod <- Lsu <- Lvar <- vector("list", B)
for (b in seq_len(B)) {
  i <- sample(length(d$time), round(0.8 * length(d$time)))
  Xs <- X[i, , drop = FALSE]; ts <- d$time[i]; ss <- d$status[i]
  Lmod[[b]] <- sel_module(ts, ss, Xs)
  Lsu[[b]]  <- survfcbf_select(ts, ss, Xs, q = 2L, top = K)$selected
  Lvar[[b]] <- cut_ranking(apply(Xs, 2L, stats::var), k = K)
}
f <- integer(length(big)); for (s in Lmod) f[s] <- f[s] + 1L
res <- data.frame(cohort = co, n = length(d$time), events = sum(d$status), p = ncol(X),
  n_modules = length(big),
  stab_module = stabm::stabilityNogueira(Lmod, p = length(big)),
  stab_gene   = stabm::stabilityNogueira(Lsu,  p = ncol(X)),
  stab_var    = stabm::stabilityNogueira(Lvar, p = ncol(X)),
  module_consensus = sum(f / B >= 0.5), stringsAsFactors = FALSE)
saveRDS(res, out)
cat(sprintf("%s: module %.3f | gene %.3f | variance %.3f | consensus %d/%d modules\n",
            co, res$stab_module, res$stab_gene, res$stab_var, res$module_consensus, K))
cat("MODULE_BOMMERT_DONE\n")
