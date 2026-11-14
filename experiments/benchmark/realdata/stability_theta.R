# Sensitivity of the stability index to its similarity threshold theta. The benchmark scores
# stability with stabm's similarity-adjusted intersection count at theta = 0.9, so only
# near-copies (|r| >= 0.9) count as interchangeable. A lower theta also credits swapping a
# gene for a correlated substitute, the representative instability that FCBF is known for.
# This rescores the ten-gene panels that panel_redundancy.R saved (<cohort>_panels.csv),
# with the same outcome-free similarity matrix and seed, at theta = 0.9 (which must
# reproduce the stability of that run) and 0.7, both with stabm's default 10,000 draws for
# the chance correction, and at theta = 0.5 with 1,000 draws. On the TCGA coding genes the
# correction at theta = 0.7 takes 5 to 30 minutes per panel and at 0.5 several hours with
# the default draws (a first run with one task per cohort timed out), so theta = 0.5 uses
# fewer draws and is the approximate one, and every (cohort, panel) is its own task.
#
# SLURM array over the 77 (cohort, panel) units, each writing its own file:
#   sbatch -p short --qos=short --time=03:00:00 --mem=32G --cpus-per-task=1 --array=1-77 \
#     cesga/submit_script.sbatch realdata/stability_theta.R
# Reads PANELS_DIR (default results/realdata_redundancy_v3), writes <cohort>_<panel>.csv to
# RD_OUTDIR (default results/realdata_stability_theta).
suppressMessages({ library(survfcbf); library(mlr3); library(mlr3proba); library(stabm)
                   library(farff); library(jsonlite); library(curl) })
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(getwd(), "realdata/Data"))
SELS <- c("survfcbf_q2_K10", "mrmr_K10", "mrmr_q2_K10", "mim_surv_K10", "mim_surv_q2_K10",
          "carss_K10", "variance_K10")
units <- expand.grid(selector = SELS, cohort = names(BOMMERT_IDS), stringsAsFactors = FALSE)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
co <- units$cohort[id]; s <- units$selector[id]
panels_dir <- Sys.getenv("PANELS_DIR", "results/realdata_redundancy_v3")
outdir <- Sys.getenv("RD_OUTDIR", "results/realdata_stability_theta")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, sprintf("%s_%s.csv", co, s))
if (file.exists(out)) { cat(basename(out), "exists, skip\n"); quit(save = "no") }
t0 <- proc.time()[["elapsed"]]

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
task <- surv_task(d$time, d$status, X)
ex <- survfcbf:::extract_surv(task)                 # the feature order the benchmark indexes
S <- abs(stats::cor(ex$X)); S[is.na(S)] <- 0; diag(S) <- 1   # as similarity_matrix(scope = "full")
rm(X, d); invisible(gc())

pn <- utils::read.csv(file.path(panels_dir, paste0(co, "_panels.csv")), stringsAsFactors = FALSE)
main <- utils::read.csv(file.path(panels_dir, paste0(co, ".csv")), stringsAsFactors = FALSE)
z <- pn[pn$selector == s, ]; z <- z[order(z$fold, z$position), ]
sets <- lapply(split(z$feature, z$fold), function(f) as.integer(match(f, ex$features)))
stopifnot("panel feature not in the task" = !anyNA(unlist(sets)))
sets <- sets[vapply(sets, length, 1L) > 0L]
GRID <- data.frame(theta = c(0.9, 0.7, 0.5), draws = c(10000L, 10000L, 1000L))
rows <- lapply(seq_len(nrow(GRID)), function(i) {
  t1 <- proc.time()[["elapsed"]]
  v <- if (length(sets) >= 2L) survfcbf:::with_seed(1L, stabm::stabilityIntersectionCount(
    sets, sim.mat = S, threshold = GRID$theta[i], N = GRID$draws[i])) else NA_real_
  cat(sprintf("%s %s theta %.1f (N %d): %.3f (%.0f s)\n", co, s, GRID$theta[i], GRID$draws[i], v,
              proc.time()[["elapsed"]] - t1))
  data.frame(cohort = co, selector = s, theta = GRID$theta[i], draws = GRID$draws[i], stability = v,
             stability_main_run = main$stability[main$selector == s], seconds = proc.time()[["elapsed"]] - t1)
})
res <- do.call(rbind, rows)
stopifnot("theta = 0.9 must reproduce the panel run" =
            isTRUE(all.equal(res$stability[1], res$stability_main_run[1], tolerance = 1e-8)))
utils::write.csv(res, out, row.names = FALSE)
cat(sprintf("%s %s done in %.1f min\nSTABILITY_THETA_DONE\n", co, s, (proc.time()[["elapsed"]] - t0) / 60))
