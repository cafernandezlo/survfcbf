# CAR-relevance ablation (the mechanism evidence). Isolates the RELEVANCE
# metric: keep the exact same machinery (correlation predominance pruning at theta=0.9 +
# top-K, on coding-only features) and swap ONLY the ranking -- SU relevance (survfcbf's
# vehicle) vs the carss IPCW CAR score. If CAR-predominance is MORE STABLE across
# subsamples than SU-predominance at matched K, the instability comes from SU, not the
# predominance pruning (the method's contribution). Per-cohort subsample-stability
# probe, SLURM array. Stability = Nogueira index (no p x p matrix needed, cheap) plus the
# selection-frequency consensus, both comparing SU vs CAR internally.
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl); library(stabm) })
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R")
Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
reg <- baseline_registry()
ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot("SLURM_ARRAY_TASK_ID out of range" = !is.na(id) && id >= 1L && id <= length(cohorts))
co <- cohorts[id]

outdir <- "results/realdata_car_ablation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".rds"))
if (file.exists(out)) { cat(co, "exists, skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
keep <- grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc
X <- X[, keep, drop = FALSE]; nm <- colnames(X)
B <- 30L; TOP <- 10L; THETA <- 0.9

## rank by score (desc), correlation-prune (|cor|>theta -> redundant, keep higher score), cap at top
pred_prune <- function(score, Xsub, top = TOP, theta = THETA) {
  ord <- order(score, decreasing = TRUE); kept <- integer(0)
  for (j in ord) {
    if (!length(kept)) { kept <- j; next }
    if (max(abs(stats::cor(Xsub[, j], Xsub[, kept, drop = FALSE]))) < theta) kept <- c(kept, j)
    if (length(kept) >= top) break
  }
  kept
}
su_relevance <- function(t, s, Xsub) {
  ym <- survival_target(t, s, method = "martingale")
  yb <- discretize_quantile(ym, q = 2L)$binned
  vapply(seq_len(ncol(Xsub)), function(j)
    symmetrical_uncertainty(discretize_quantile(Xsub[, j], q = 2L)$binned, yb), numeric(1))
}
car_relevance <- function(t, s, Xsub) run_baseline(ent("carss"), t, s, Xsub, prop = NULL, seed = 1L)$score

set.seed(1)
sel_su <- sel_car <- vector("list", B)
for (b in seq_len(B)) {
  idx <- sample(length(d$time), round(0.8 * length(d$time)))
  Xs <- X[idx, , drop = FALSE]; ts <- d$time[idx]; ss <- d$status[idx]
  sel_su[[b]]  <- pred_prune(su_relevance(ts, ss, Xs),  Xs)
  sel_car[[b]] <- pred_prune(car_relevance(ts, ss, Xs), Xs)
}
p <- ncol(X)
freq <- function(sel) { f <- integer(p); for (s in sel) f[s] <- f[s] + 1L; f / B }
fsu <- freq(sel_su); fcar <- freq(sel_car)
res <- data.frame(
  cohort = co, n = length(d$time), events = sum(d$status), p = p,
  stab_su  = stabm::stabilityNogueira(sel_su,  p = p),
  stab_car = stabm::stabilityNogueira(sel_car, p = p),
  maxfreq_su = max(fsu), maxfreq_car = max(fcar),
  consensus_su = sum(fsu >= 0.5), consensus_car = sum(fcar >= 0.5),
  stringsAsFactors = FALSE)
saveRDS(res, out)
cat(sprintf("%s: stab SU %.3f vs CAR %.3f | maxfreq SU %.2f vs CAR %.2f | consensus SU %d vs CAR %d\n",
            co, res$stab_su, res$stab_car, res$maxfreq_su, res$maxfreq_car, res$consensus_su, res$consensus_car))
cat("CAR_ABLATION_DONE\n")
