# Combine the per-cohort CSVs from the real-data array into one
# table. Run from experiments/benchmark inside the container:
#   apptainer exec "$SIF" Rscript --no-init-file <abs>/realdata/aggregate_realdata.R
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench)) setwd(bench)

cells_dir <- "results/realdata/cells"
files <- list.files(cells_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("no cohort CSVs in ", cells_dir, call. = FALSE)

raw <- do.call(rbind, lapply(files, utils::read.csv, check.names = FALSE))
out <- file.path("results", "benchmark_realdata.csv")
utils::write.csv(raw, out, row.names = FALSE)
cat(sprintf("%d / 11 cohorts present -> %s (%d selector-rows)\n", length(files), out, nrow(raw)))
print(table(raw$dataset))
