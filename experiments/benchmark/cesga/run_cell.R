# CESGA SLURM array driver. Each array task processes a STRIDED
# subset of grid cells -- cell indices i with (i-1) %% N_TASKS == (TASK_ID-1) --
# so heavy and cheap scenarios spread evenly across tasks and the array stays
# under the QoS MaxSubmitPU cap (short allows 100 queued jobs/user, a 330-cell
# 1-task-per-cell array is rejected). With N_TASKS unset each task does one cell.
# Cells already written are skipped, so a timed-out task can simply be re-run.
# Run inside the bioconductor R 4.5.2 apptainer container (see submit_grid.sbatch).
# Args: <tier>.
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners); library(stabm)
})
bench <- Sys.getenv("BENCH_DIR", unset = "")     # set the working dir explicitly (apptainer cwd is unreliable)
if (nzchar(bench)) setwd(bench)
source("grid_lib.R")

tier <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(tier)) tier <- "cesga_core"
cfg <- grid_config(tier)
scn <- build_scenarios(cfg)
cells <- grid_cells(cfg)
n_cells <- nrow(cells)

id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_tasks <- as.integer(Sys.getenv("N_TASKS", as.character(n_cells)))
stopifnot("N_TASKS must be a positive integer" = !is.na(n_tasks) && n_tasks >= 1L,
          "SLURM_ARRAY_TASK_ID out of range" = !is.na(id) && id >= 1L && id <= n_tasks)

mine <- seq.int(id, n_cells, by = n_tasks)   # strided assignment: cells id, id+n_tasks, ...
outdir <- file.path("results", tier, "cells")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (k in mine) {
  out <- file.path(outdir, sprintf("cell_%05d.csv", k))
  if (file.exists(out)) { cat(sprintf("cell %d/%d skipped (exists)\n", k, n_cells)); next }
  g <- cells[k, ]
  t0 <- proc.time()[["elapsed"]]
  res <- run_cell(scn[[g$s]], g$learner, g$rep, cfg, master_seed = 1L)
  utils::write.csv(res, out, row.names = FALSE)
  cat(sprintf("cell %d/%d done: %s / %s / rep %d -> %d rows, %.1f min\n",
              k, n_cells, g$scenario, g$learner, g$rep, nrow(res),
              (proc.time()[["elapsed"]] - t0) / 60))
}
cat(sprintf("task %d/%d done: %d cells\n", id, n_tasks, length(mine)))
