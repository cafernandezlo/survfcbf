# Figure of the Application Note, drawn from the tracked result tables so its numbers cannot
# drift from the text. Conventions: Okabe-Ito colours assigned by entity and never cycled,
# redundant symbols and line types so the figure survives greyscale, recessive grid.
# A: what crosses the fold boundary, as nested_cv_select() and score_benchmark() do it.
# B: panel size per cohort (a dot strip, since bars on a log axis have no meaningful baseline).
# C: selection runtime against p (synthetic data, realdata/runtime_scaling.R).
#
# Usage (from the repository root):  Rscript --no-init-file paper/make_appnote_figure.R
# Output: paper/appnote/Figure_appnote.png at 300 dpi, double-column width.
root <- Sys.getenv("SURVFCBF_ROOT", unset = "."); setwd(root)
RES <- "experiments/benchmark/results"

OI <- c(blue = "#0072B2", green = "#009E73", orange = "#E69F00", vermillion = "#D55E00",
        purple = "#CC79A7")
INK <- "grey25"; RANKERS <- "grey35"
big <- function(x) formatC(x, format = "d", big.mark = ",")

grDevices::png("paper/appnote/Figure_appnote.png", width = 180, height = 118, units = "mm",
               res = 300, type = "cairo", pointsize = 9)
graphics::layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = TRUE), heights = c(0.82, 1))
AXES <- list(mgp = c(2.2, 0.7, 0), tcl = -0.25, cex.axis = 0.85, cex.lab = 0.95, bty = "n", las = 1)
heading <- function(letter, s) graphics::mtext(paste0(letter, "   ", s), side = 3, adj = 0, line = 0.6,
                                               font = 2, cex = 0.85)

## ---------------------------------------------------------------- A: the fold boundary
do.call(graphics::par, c(list(mar = c(0.2, 2.2, 2.2, 0.8)), AXES))
graphics::plot.new(); graphics::plot.window(xlim = c(0, 100), ylim = c(0, 40), xaxs = "i", yaxs = "i")
heading("A", "Everything learned from data stays inside the training fold")

graphics::rect(0, 21.5, 100, 40, col = "grey95", border = NA)
graphics::rect(0, 0, 100, 16.5, col = "white", border = "grey80", lwd = 0.6)
graphics::text(-1.1, 30.7, "Training", srt = 90, cex = 0.8, font = 2, col = INK, xpd = NA)
graphics::text(-1.1, 8.2, "Held-out", srt = 90, cex = 0.8, font = 2, col = INK, xpd = NA)

step <- function(x, y, label, w = 14.6, h = 9.6, border = "grey55", lwd = 0.7) {
  graphics::rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = "white", border = border, lwd = lwd)
  graphics::text(x, y, label, cex = 0.76, col = "black")
}
flow <- function(x0, y0, x1, y1, lwd = 0.9, col = "grey35")
  graphics::arrows(x0, y0, x1, y1, length = 0.05, lwd = lwd, col = col)

xs <- 9.3 + 16.3 * 0:5; yt <- 28.6
labels <- c("Null Cox model:\nmartingale residuals", "Equal-frequency\nbin cut-points",
            "SU relevance\nwith the residuals", "Markov-blanket\npruning", "Panel of\nK genes",
            "Ridge Cox\nfitted on the panel")
for (i in 1:6) step(xs[i], yt, labels[i], border = if (i %in% 3:4) OI[["vermillion"]] else "grey55",
                    lwd = if (i %in% 3:4) 1.1 else 0.7)
for (i in 1:5) flow(xs[i] + 7.3, yt, xs[i + 1] - 7.3, yt)
graphics::segments(xs[1] - 7.3, 35.6, xs[5] + 7.3, 35.6, col = "grey55", lwd = 0.7)
graphics::text(mean(xs[c(1, 5)]), 37.6, "survfcbf_select(), on the training rows only", cex = 0.76, col = INK)
graphics::text(xs[6], 37.6, "any learner", cex = 0.76, col = INK)

## the boundary: the only things allowed across, all in one direction
graphics::segments(0, 19, 100, 19, lty = 2, lwd = 0.9, col = "grey25")
graphics::text(1, 19, "fold boundary", adj = c(0, -0.35), cex = 0.72, col = INK, font = 3)
yb <- 8.2
step(xs[6] - 8.15, yb, "Predict on the held-out samples,\nusing the panel genes only", w = 30.9, h = 9)
step(xs[4] - 8.15, yb, "Score: truncated Uno's C,\nintegrated Brier score", w = 30.9, h = 9)
flow(xs[4] - 8.15, 21.5, xs[4] - 8.15, yb + 4.5, lwd = 1.3, col = "black")
flow(xs[5], yt - 4.8, xs[5], yb + 4.5, lwd = 1.3, col = "black")
flow(xs[6], yt - 4.8, xs[6], yb + 4.5, lwd = 1.3, col = "black")
flow(xs[6] - 8.15 - 15.45, yb, xs[4] - 8.15 + 15.45, yb)
graphics::text(xs[4] - 8.95, 19, "censoring weights, truncation time", adj = c(1, 1.25), cex = 0.72, col = "black")
graphics::text(xs[5] - 0.8, 19, "gene names", adj = c(1, 1.25), cex = 0.72, col = "black")
graphics::text(xs[6] - 0.8, 19, "fitted model", adj = c(1, 1.25), cex = 0.72, col = "black")
graphics::text(xs[1] + 2.5, yb, "Repeated for each outer fold.\nHeld-out samples never reach\nthe steps above the line.",
               cex = 0.74, col = INK)

## ---------------------------------------------------------------- B: panel size
rd <- utils::read.csv(file.path(RES, "benchmark_realdata.csv"), check.names = FALSE)
pars <- function(s) rd$parsimony[rd$selector == s & !is.na(rd$parsimony)]
## the three rankers keep the same 10% of the features, so they share one row
stopifnot(identical(pars("variance"), pars("carss")), identical(pars("variance"), pars("mim_surv")))
rows <- list(
  list(label = "Rankers at 10%\n(variance, carss, MI)", v = pars("variance"), col = RANKERS, pch = 17, f = big),
  list(label = "Penalised Cox\n(lasso)", v = pars("cox_glmnet_l1"), col = OI[["purple"]], pch = 15, f = function(x) formatC(x, format = "f", digits = 1)),
  list(label = "survfcbf,\ntop-K mode", v = pars("survfcbf_topk"), col = OI[["vermillion"]], pch = 19, f = function(x) formatC(x, format = "d")),
  list(label = "survfcbf,\npermutation mode", v = pars("survfcbf"), col = OI[["vermillion"]], pch = 1, f = function(x) formatC(x, format = "f", digits = 1)))
do.call(graphics::par, c(list(mar = c(3.6, 8.2, 2.2, 1.6)), AXES))
ys <- rev(seq_along(rows))
graphics::plot(NA, xlim = c(0.8, 9000), ylim = c(0.5, length(rows) + 0.5), log = "x", axes = FALSE,
               xlab = "features selected per cohort (log scale)", ylab = "")
graphics::abline(v = 10^(0:3), col = "grey88", lwd = 0.6)
graphics::axis(1, at = 10^(0:3), labels = c("1", "10", "100", "1,000"), lwd = 0.6, col = "grey55")
graphics::axis(2, at = ys, labels = vapply(rows, `[[`, "", "label"), lwd = 0, cex.axis = 0.82, line = -0.6)
for (i in seq_along(rows)) { r <- rows[[i]]
  set.seed(i); jit <- stats::runif(length(r$v), -0.13, 0.13)
  graphics::points(r$v, ys[i] + jit, pch = if (r$pch == 1) 1 else 16, cex = 0.55,
                   col = grDevices::adjustcolor(r$col, alpha.f = 0.55))
  m <- mean(r$v)
  graphics::segments(m, ys[i] - 0.3, m, ys[i] + 0.3, col = r$col, lwd = 2)
  graphics::text(m, ys[i] + 0.3, r$f(m), pos = 3, offset = 0.2, cex = 0.82, col = "black", xpd = NA)
}
heading("B", "Panel size on the eleven TCGA cohorts")

## ---------------------------------------------------------------- C: runtime
rt <- utils::read.csv(file.path(RES, "realdata_runtime/runtime_scaling.csv"), check.names = FALSE)
slope <- unname(stats::coef(stats::lm(log(survfcbf_topk) ~ log(p), data = rt))[2])
secs  <- rt$survfcbf_topk[which.max(rt$p)]
series <- list(survfcbf_topk = list("survfcbf top-K", OI[["vermillion"]], 19, 1),
               carss    = list("carss",    OI[["green"]],  15, 3),
               mim_surv = list("MI ranker", OI[["orange"]], 18, 4),
               variance = list("variance", OI[["blue"]],   17, 2))
do.call(graphics::par, c(list(mar = c(3.6, 3.8, 2.2, 1.8)), AXES))
graphics::plot(NA, xlim = range(rt$p), ylim = c(0.008, 60), log = "xy", axes = FALSE,
               xlab = "features p (log scale)", ylab = "selection time, s (log scale)")
graphics::abline(h = 10^(-2:1), v = c(1000, 2000, 5000, 10000, 20000), col = "grey88", lwd = 0.6)
graphics::axis(1, at = c(1000, 2000, 5000, 10000, 20000), labels = big(c(1000, 2000, 5000, 10000, 20000)),
               lwd = 0.6, col = "grey55")
graphics::axis(2, at = 10^(-2:1), labels = c("0.01", "0.1", "1", "10"), lwd = 0.6, col = "grey55")
for (nm in names(series)) { z <- series[[nm]]
  graphics::lines(rt$p, rt[[nm]], type = "b", col = z[[2]], pch = z[[3]], lty = z[[4]], lwd = 1.6, cex = 0.8) }
graphics::text(max(rt$p), secs * 1.35, bquote(p^.(formatC(slope, format = "f", digits = 2)) * ", " *
               .(formatC(secs, format = "f", digits = 1)) * " s"), adj = c(1, 0), cex = 0.82, col = "black")
graphics::legend("topleft", vapply(series, `[[`, "", 1), col = vapply(series, `[[`, "", 2),
                 pch = vapply(series, `[[`, 0, 3), lty = vapply(series, `[[`, 0, 4),
                 lwd = 1.6, bty = "n", cex = 0.78, bg = "white", y.intersp = 0.95)
heading("C", "Selection time is close to linear in p")
grDevices::dev.off()
cat("wrote paper/appnote/Figure_appnote.png\nAPPNOTE_FIGURE_DONE\n")
