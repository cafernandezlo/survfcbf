# Simulation benchmark grid -- local targets orchestration.
#
# The grid definition (scenarios, selectors, cell runner, aggregation) lives in
# grid_lib.R, shared with the CESGA SLURM array (cesga/run_cell.R) so there is a
# single source of truth. `.tier` picks the configuration (grid_config()):
# "demo" validates locally; "cesga_core" / "production" are the HPC tiers (run via
# the array, not this serial pipeline). Run with targets::tar_make().

library(targets)
tar_option_set(packages = c("survfcbf", "mlr3", "mlr3proba", "mlr3extralearners", "stabm"))
source("grid_lib.R")

.tier <- "demo"
cfg <- grid_config(.tier)

build_raw <- function(cfg, master_seed) {
  scn <- build_scenarios(cfg)
  cells <- grid_cells(cfg)
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    g <- cells[i, ]
    run_cell(scn[[g$s]], g$learner, g$rep, cfg, master_seed)
  }))
}

list(
  tar_target(master_seed, 1L),
  tar_target(raw_results, build_raw(cfg, master_seed)),
  tar_target(results_table, aggregate_raw(raw_results)),
  tar_target(results_csv, {
    dir.create("results", showWarnings = FALSE)
    path <- file.path("results", sprintf("benchmark_%s.csv", .tier))
    utils::write.csv(results_table, path, row.names = FALSE)
    path
  }, format = "file")
)
