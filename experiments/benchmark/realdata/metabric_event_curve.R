# METABRIC event-downsampling stability curve (the real-data test of the
# events-dependency mechanism found in simulation). METABRIC has 1143
# events; we subsample it down to fewer events (keeping the natural censoring ratio, so a
# 75-event draw is a ~130-patient PAAD-like cohort) and measure the subsample stability of
# each selector's top-10 at every event level. The simulation predicts survfcbf's SU
# instability worsens as events drop below ~300-400; if that curve reproduces on real
# METABRIC (and variance stays flat-high), the mechanism is confirmed on real data.
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench) && dir.exists(bench)) setwd(bench)
## load survfcbf from source (the installed pkg may predate the top-K `top` arg)
if (requireNamespace("survfcbf", quietly = TRUE) && "top" %in% names(formals(survfcbf::survfcbf_select))) {
  suppressMessages(library(survfcbf))
} else {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))
}
if (requireNamespace("stabm", quietly = TRUE)) library(stabm)
source("realdata/load_metabric.R")
set.seed(1)

d <- load_metabric()
cat(sprintf("METABRIC loaded: n=%d events=%d p=%d (log2 microarray)\n",
            length(d$time), sum(d$status), ncol(d$X)))
Es  <- c(75L, 150L, 300L, 500L, 800L, sum(d$status))
R   <- 3L; B <- 20L; TOP <- 10L
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
selfun <- list(
  survfcbf_topk = function(t, s, X) survfcbf_select(t, s, X, q = 2L, top = TOP)$selected,
  variance      = function(t, s, X) cut_ranking(run_baseline(ent("variance"), t, s, X, prop = NULL, seed = 1L)$score, k = TOP),
  carss         = function(t, s, X) cut_ranking(run_baseline(ent("carss"),    t, s, X, prop = NULL, seed = 1L)$score, k = TOP)
)
## subsample to ~E events, keeping the natural censoring ratio (fewer events -> smaller cohort)
sub_to_events <- function(status, E, seed) {
  set.seed(seed); ev <- which(status == 1L); cn <- which(status == 0L)
  ne <- min(E, length(ev)); nc <- min(length(cn), round(length(cn) * ne / length(ev)))
  c(sample(ev, ne), sample(cn, nc))
}

rows <- list()
for (E in Es) for (r in seq_len(R)) {
  idx <- sub_to_events(d$status, E, seed = E * 100L + r)
  ts <- d$time[idx]; ss <- d$status[idx]; Xs <- d$X[idx, , drop = FALSE]
  for (sl in names(selfun)) {
    sel <- vector("list", B)
    for (b in seq_len(B)) {
      set.seed(E * 1000L + r * 10L + b)
      ii <- sample(length(ts), round(0.8 * length(ts)))
      sel[[b]] <- tryCatch(selfun[[sl]](ts[ii], ss[ii], Xs[ii, , drop = FALSE]),
                           error = function(e) integer(0))
    }
    rows[[length(rows) + 1L]] <- data.frame(events = E, n = length(idx), rep = r, selector = sl,
      stability = stabm::stabilityNogueira(sel, p = ncol(Xs)), stringsAsFactors = FALSE)
  }
  cat(sprintf("  E=%4d (n=%4d) rep %d done\n", E, length(idx), r))
}
res <- do.call(rbind, rows)
OUT <- "results/realdata_metabric"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(OUT, "event_curve_raw.csv"), row.names = FALSE)
agg <- aggregate(stability ~ events + selector, res, function(x) mean(x, na.rm = TRUE))
wide <- reshape(agg, idvar = "events", timevar = "selector", direction = "wide")
utils::write.csv(wide, file.path(OUT, "event_curve.csv"), row.names = FALSE)
cat("\n=== stability vs events (mean over reps) ===\n"); print(wide, row.names = FALSE)

grDevices::png(file.path(OUT, "event_curve.png"), width = 900, height = 620, res = 130)
op <- graphics::par(mar = c(4.5, 4.5, 3, 1))
cols <- c(survfcbf_topk = "firebrick", variance = "grey40", carss = "steelblue")
sels <- names(cols)
graphics::plot(NA, xlim = range(Es), ylim = c(0, 1), log = "x", xlab = "events (subsampled)",
               ylab = "selection stability (Nogueira)", main = "METABRIC: stability vs events (top-10)")
for (sl in sels) { a <- agg[agg$selector == sl, ]; a <- a[order(a$events), ]
  graphics::lines(a$events, a$stability, type = "b", pch = 19, col = cols[sl]) }
graphics::abline(v = 350, lty = 3, col = "grey70")
graphics::legend("topleft", c(sels, "sim crossover ~350 ev"), col = c(cols, "grey70"),
                 pch = c(19, 19, 19, NA), lty = c(1, 1, 1, 3), bty = "n")
graphics::par(op); grDevices::dev.off()
cat(sprintf("\ntables + figure -> %s\nMETABRIC_CURVE_DONE\n", OUT))
