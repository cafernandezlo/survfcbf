# Combine the per-cell CSVs from the SLURM array into the raw and aggregated
# benchmark tables. Run from experiments/benchmark inside the
# container: apptainer exec "$SIF" Rscript --no-init-file <abs>/cesga/aggregate.R <tier>
suppressMessages(library(survfcbf))              # grid_lib.R's grid_cells() needs default_params()
bench <- Sys.getenv("BENCH_DIR", unset = "")     # explicit working dir (apptainer cwd is unreliable)
if (nzchar(bench)) setwd(bench)
source("grid_lib.R")

tier <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(tier)) tier <- "cesga_core"

cells_dir <- file.path("results", tier, "cells")
files <- list.files(cells_dir, pattern = "^cell_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("no cell CSVs found in ", cells_dir, call. = FALSE)

raw <- do.call(rbind, lapply(files, utils::read.csv, check.names = FALSE))
expected <- nrow(grid_cells(grid_config(tier)))
cat(sprintf("%d / %d cells present (%d selector-rows)\n", length(files), expected, nrow(raw)))

utils::write.csv(raw, file.path("results", sprintf("raw_%s.csv", tier)), row.names = FALSE)
agg <- aggregate_raw(raw)
utils::write.csv(agg, file.path("results", sprintf("benchmark_%s.csv", tier)), row.names = FALSE)
cat(sprintf("wrote raw_%s.csv (%d rows) and benchmark_%s.csv (%d rows)\n",
            tier, nrow(raw), tier, nrow(agg)))
