# Uno's C truncated at a fold's own follow-up horizon. Without a truncation time, an event
# in the test fold later than the last censoring of the training fold meets a censoring
# survival estimate of zero, its pairs take the capped weight 1 / eps^2, and a handful of
# pairs can carry nearly all the weight of the fold (found on METABRIC, where one fold
# scores about 0.97 for every model). Truncating at tau = the p-quantile of the training
# fold's observed times keeps the estimate on the range where the censoring distribution
# is supported (Uno et al. 2011). tau is learned on the training rows only.
#
# uno_by_fold(bm) returns, per scored fold, tau, the untruncated Uno's C and the truncated
# one, so a degenerate fold can be seen rather than averaged away.
uno_by_fold <- function(bm, p = 0.8) {
  do.call(rbind, lapply(seq_along(bm$folds), function(f) {
    fd <- bm$folds[[f]]
    if (is.null(fd$prediction)) return(NULL)
    tt <- bm$task$truth(rows = fd$train_set)[, 1]
    tau <- unname(stats::quantile(tt, p))
    m0 <- mlr3::msr("surv.cindex", weight_meth = "G2", id = "uno_c")
    m1 <- mlr3::msr("surv.cindex", weight_meth = "G2", t_max = tau, id = "uno_c_trunc")
    v <- fd$prediction$score(list(m0, m1), task = bm$task, train_set = fd$train_set)
    data.frame(fold = f, tau = tau, uno_c = v[["uno_c"]], uno_c_trunc = v[["uno_c_trunc"]])
  }))
}
