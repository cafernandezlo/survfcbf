# Lightweight CESGA env smoke test: confirms the installed
# survfcbf + mlr3 stack run benchmark_selectors end-to-end inside the container,
# before submitting the full array. Run with the absolute path:
#   apptainer exec "$SIF" Rscript --no-init-file <abs>/cesga/cesga_smoke.R
suppressMessages({
  library(survfcbf); library(mlr3); library(mlr3proba); library(mlr3extralearners); library(stabm)
})
if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")

sim <- simulate_dataset(default_params(p = 40, n = 150, seed = 1, calib_n_presample = 400,
        k_drivers = 2, cluster_size = 3, n_distractor_subblocks = 2, distractor_size = 3, target_C = 0.90))
task <- surv_task(sim$time, sim$status, sim$X)
sels <- list(variance = function(t, s, X) variance_select(t, s, X, prop = 0.1)$selected,
             survfcbf = function(t, s, X) survfcbf_select(t, s, X, q = 8, B = 20)$selected)
tab <- benchmark_selectors(task, sels, mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min"),
        outer_folds = 3, manifest = sim$manifest, seed = 1)
print(tab[, c("selector", "uno_c", "stability", "parsimony", "sensitivity")], digits = 3)
cat("SMOKE_OK\n")
