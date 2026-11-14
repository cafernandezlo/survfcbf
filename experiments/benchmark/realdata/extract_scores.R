# Extract the FULL per-gene relevance score of every selector for a GSEA analysis
# (interpretability pillar, the ranking-based complement to the top-K
# ORA). GSEA tests the whole ranking, not an arbitrary top-K cut, and so asks the deeper
# mechanistic question: does each selector's relevance metric concentrate cancer/hallmark
# biology at the top? survfcbf is ranked by its SU relevance (SU of each gene with the
# martingale target at q=2, the top-K config); the baselines by their own per-gene score.
# Coding-only (protein-coding rna features), so genes map cleanly to the MSigDB symbol
# sets. SLURM array, one cohort/task. Run in the apptainer container on CESGA.
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl) })
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

outdir <- Sys.getenv("RD_SCORES_DIR", "results/realdata_scores")   # a rerun writes elsewhere
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".rds"))
if (file.exists(out)) { cat(co, "exists, skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
keep <- grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc      # coding rna only
X <- X[, keep, drop = FALSE]; nm <- colnames(X)

## survfcbf relevance = SU(gene, martingale target) at q=2 (the top-K config)
ym <- survival_target(d$time, d$status, method = "martingale")
yb <- discretize_quantile(ym, q = 2L)$binned
su <- vapply(seq_len(ncol(X)), function(j)
  symmetrical_uncertainty(discretize_quantile(X[, j], q = 2L)$binned, yb), numeric(1))

scores <- data.frame(
  feature  = nm,
  survfcbf = su,
  ## the secondary key survfcbf_select() uses to order genes tied on SU
  survfcbf_tiebreak = survfcbf:::spearman_tiebreak(X, ym),
  variance = run_baseline(ent("variance"), d$time, d$status, X, prop = NULL, seed = 1L)$score,
  carss    = run_baseline(ent("carss"),    d$time, d$status, X, prop = NULL, seed = 1L)$score,
  mim_surv = run_baseline(ent("mim_surv"), d$time, d$status, X, prop = NULL, seed = 1L)$score,
  stringsAsFactors = FALSE)
saveRDS(scores, out)

## What the q = 2 target holds: the upper half of the martingale residuals takes every
## event and, when events are fewer than half the patients, the censored patients with
## the shortest follow-up. Recorded so the supplement can report it per cohort.
upper <- yb == max(yb)
cens <- d$status == 0
target <- data.frame(
  cohort = co, n = length(d$time), events = sum(d$status), upper_n = sum(upper),
  upper_events = sum(upper & !cens), upper_censored = sum(upper & cens),
  followup_censored_upper = stats::median(d$time[upper & cens]),
  followup_censored_lower = stats::median(d$time[!upper & cens]))
tdir <- file.path(dirname(outdir), "realdata_target")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(target, file.path(tdir, paste0(co, ".csv")), row.names = FALSE)
cat(sprintf("%s scores done: %d coding features\n", co, nrow(scores)))
cat("SCORES_DONE\n")
