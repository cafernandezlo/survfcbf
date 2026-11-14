# Runtime scaling (the honest efficiency read-out). Times the selection
# step of survfcbf_topk and the filter baselines as p grows (fixed n), to characterise
# FCBF's classic near-linearity claim -- and to report it honestly (top-K is far faster
# than its own delta / cv_glmnet, but the trivial filters are faster still). Fits the
# scaling exponent runtime ~ p^alpha for survfcbf_topk. Local, synthetic data.
setwd(Sys.getenv("BENCH_DIR", unset = "."))
suppressMessages(devtools::load_all("../..", quiet = TRUE))   # from experiments/benchmark
set.seed(1)

ps <- c(1000L, 2000L, 5000L, 10000L, 20000L)
n  <- 300L
## median of REPS wall-clock runs, so one busy moment on the machine cannot distort a point
REPS <- 5L
timeit <- function(expr) {
  f <- substitute(expr); env <- parent.frame()
  stats::median(vapply(seq_len(REPS), function(i) {
    t0 <- proc.time()[["elapsed"]]; eval(f, env); proc.time()[["elapsed"]] - t0 }, numeric(1)))
}
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]

rows <- list()
for (p in ps) {
  sim <- simulate_dataset(default_params(p = p, n = n, k_drivers = 3, cluster_size = 5,
          n_distractor_subblocks = 4, distractor_size = 8, censor_fraction = 0.7,
          target_C = 0.7, seed = 1, calib_n_presample = 800))
  ti <- sim$time; st <- sim$status; X <- sim$X
  t_topk <- timeit(survfcbf_select(ti, st, X, q = 2L, top = 10L))
  t_var  <- timeit(run_baseline(ent("variance"), ti, st, X, prop = 0.1, seed = 1L))
  t_car  <- timeit(run_baseline(ent("carss"),    ti, st, X, prop = 0.1, seed = 1L))
  t_mim  <- timeit(run_baseline(ent("mim_surv"), ti, st, X, prop = 0.1, seed = 1L))
  rows[[length(rows) + 1L]] <- data.frame(p = p, survfcbf_topk = t_topk, variance = t_var,
                                          carss = t_car, mim_surv = t_mim)
  cat(sprintf("p=%6d | topk %.2fs  var %.2fs  carss %.2fs  mim %.2fs\n", p, t_topk, t_var, t_car, t_mim))
}
tab <- do.call(rbind, rows)
OUT <- "results/realdata_runtime"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(tab, file.path(OUT, "runtime_scaling.csv"), row.names = FALSE)

## scaling exponent for survfcbf_topk (log-log slope)
fit <- stats::lm(log(survfcbf_topk) ~ log(p), data = tab)
alpha <- unname(coef(fit)[2])
cat(sprintf("\nsurvfcbf_topk scaling: runtime ~ p^%.2f (alpha~1 = near-linear)\n", alpha))

grDevices::png(file.path(OUT, "runtime_scaling.png"), width = 900, height = 620, res = 130)
op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
cols <- c(survfcbf_topk = "firebrick", variance = "grey50", carss = "steelblue", mim_surv = "darkgreen")
graphics::plot(tab$p, tab$survfcbf_topk, log = "xy", type = "b", pch = 19, col = cols["survfcbf_topk"],
               xlab = "p (features)", ylab = "selection time (s)", ylim = range(unlist(tab[-1])),
               main = sprintf("Runtime vs p (n=%d)  |  survfcbf_topk ~ p^%.2f", n, alpha))
for (s in c("variance", "carss", "mim_surv"))
  graphics::lines(tab$p, tab[[s]], type = "b", pch = 19, col = cols[s])
graphics::legend("topleft", names(cols), col = cols, pch = 19, lty = 1, bty = "n")
graphics::par(op); grDevices::dev.off()
cat(sprintf("tables + figure -> %s\nRUNTIME_DONE\n", OUT))
