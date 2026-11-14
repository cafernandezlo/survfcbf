# Realised-diagnostics sub-pipeline for the simulator.
#
# Runs the default scenario, the censoring sweep, and the low-variance-signal
# scenario through run_replicates() and assembles the realised-diagnostics table
# (censoring fraction, concordance, realised intra-cluster correlation, effective
# events) with Monte Carlo standard errors. Also sanity-checks the
# low-variance-signal scenario. Depends on the installed survfcbf package.
#
# Run from this directory with: targets::tar_make()

library(targets)
tar_option_set(packages = "survfcbf")

# Demonstration settings; the ficha production value is R = 200.
.R <- 50L
.master_seed <- 1L
.p <- 5000L

build_diagnostics_table <- function(R, master_seed, p) {
  scenarios <- list(
    default       = default_params(p = p, censor_fraction = 0.5),
    censoring_0.2 = default_params(p = p, censor_fraction = 0.2),
    censoring_0.7 = default_params(p = p, censor_fraction = 0.7),
    low_variance  = default_params(p = p, censor_fraction = 0.5, driver_var_scale = 0.2)
  )
  do.call(rbind, lapply(names(scenarios), function(nm) {
    s <- run_replicates(scenarios[[nm]], R = R, master_seed = master_seed)$summary
    data.frame(scenario = nm, s, row.names = NULL)
  }))
}

low_variance_check <- function(p) {
  sim <- simulate_dataset(default_params(p = p, driver_var_scale = 0.2, seed = 1L))
  m <- sim$manifest
  sig <- mean(apply(sim$X[, m$role %in% c("driver", "redundant")], 2L, stats::var))
  noise <- mean(apply(sim$X[, m$role == "noise"], 2L, stats::var))
  data.frame(signal_var = sig, noise_var = noise, ratio = sig / noise)
}

list(
  tar_target(R, .R),
  tar_target(master_seed, .master_seed),
  tar_target(p, .p),
  tar_target(diagnostics_table, build_diagnostics_table(R, master_seed, p)),
  tar_target(low_variance_sanity, low_variance_check(p)),
  tar_target(diagnostics_csv, {
    dir.create("results", showWarnings = FALSE)
    path <- file.path("results", "diagnostics_table.csv")
    utils::write.csv(diagnostics_table, path, row.names = FALSE)
    path
  }, format = "file")
)
