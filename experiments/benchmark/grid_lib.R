# Shared grid definition for the simulation benchmark.
# Pure functions, no targets dependency, sourced by both _targets.R (local
# orchestration) and cesga/run_cell.R (the SLURM array driver), so the scenario
# grid and selector set have a single source of truth.

grid_config <- function(tier = "demo") {
  cfgs <- list(
    demo = list(
      R = 1L, outer = 3L, inner_folds = 2L, B = 20L,
      p_ref = 150L, n_ref = 400L, censor_ref = 0.30,
      n_levels = c(300L, 600L), censor_levels = numeric(0), p_levels = integer(0), oat = character(0),
      q_grid = c(5L, 10L), q_fixed = c(5L, 10L), prop = 0.1, learners = "l2cox",
      selectors = c("survfcbf", "survfcbf_cvq", "variance", "cox_glmnet_l1")
    ),
    # trimmed first CESGA run: maps A + B, L2-Cox, core selectors (no grplasso /
    # at-prop / RSF-selector yet), modest R -- gets the headline maps flowing.
    # Reduced CV/null rigor (outer 5, inner 3, B 100) so survfcbf_cvq's nested CV
    # stays well under the 6h task limit on the p=2000/5000 cells; the production
    # tier below keeps the full outer=10 / inner=5 / B=200.
    cesga_core = list(
      R = 30L, outer = 5L, inner_folds = 3L, B = 100L,
      p_ref = 2000L, n_ref = 600L, censor_ref = 0.50,
      n_levels = c(150L, 300L, 600L, 1200L, 2400L), censor_levels = c(0.20, 0.50, 0.70),
      p_levels = c(200L, 1000L, 5000L), oat = character(0),
      q_grid = c(4L, 7L, 10L), q_fixed = c(4L, 7L, 10L), prop = 0.1, learners = "l2cox",
      selectors = c("survfcbf", "survfcbf_cvq", "variance", "carss", "cox_glmnet_l1", "mim_surv")
    ),
    # the full grid: all selectors, both learners, OAT axes, R=100. Dimensions are
    # genomic-realistic, matched to Bommert's 11 TCGA datasets (p=22k-27k, n=124-735,
    # events=52-152): p_ref/p_levels reach ~20k, n 125-735, censoring 0.5-0.9 lands
    # events at ~50-150. OPEN before launch: (i) cost -- cells at
    # p=20k run ~10-30h, so this needs pruning (cvq/grplasso on a cell subset, lower R
    # on secondary axes, QoS long); (ii) prefilter -- run on full p or a content-neutral
    # prefilter, NOT a variance prefilter (variance is a baseline, it would contaminate).
    production = list(
      R = 100L, outer = 10L, inner_folds = 5L, B = 200L,
      p_ref = 10000L, n_ref = 300L, censor_ref = 0.70,
      n_levels = c(125L, 250L, 500L, 735L), censor_levels = c(0.50, 0.70, 0.90),
      p_levels = c(5000L, 10000L, 20000L), oat = c("indep", "aft"),
      q_grid = c(4L, 7L, 10L), q_fixed = c(4L, 7L, 10L), prop = 0.1,
      learners = c("l2cox", "rsf"), selectors = NULL
    )
  )
  cfgs[[tier]]
}

build_scenarios <- function(cfg) {
  base <- list(k_drivers = 3L, cluster_size = 5L, n_distractor_subblocks = 4L, distractor_size = 8L,
               target_C = 0.70, p = cfg$p_ref, n = cfg$n_ref, censor_fraction = cfg$censor_ref)
  mk <- function(name, map, over) list(name = name, map = map,
          params = do.call(default_params, utils::modifyList(base, over)))
  scn <- list()
  for (n in cfg$n_levels)       scn <- c(scn, list(mk(sprintf("A_n%d", n), "A_events", list(n = n))))
  for (cs in cfg$censor_levels) scn <- c(scn, list(mk(sprintf("A_cens%02.0f", cs * 100), "A_censor", list(censor_fraction = cs))))
  for (p in cfg$p_levels)       scn <- c(scn, list(mk(sprintf("B_p%d", p), "B_dim", list(p = p))))
  if ("indep" %in% cfg$oat) scn <- c(scn, list(mk("OAT_indep", "OAT", list(cluster_size = 1L, k_drivers = 5L))))
  if ("aft" %in% cfg$oat)   scn <- c(scn, list(mk("OAT_aft", "OAT", list(time_model = "loglogistic_aft"))))
  scn
}

# l2cox uses lambda.min: the default 1se rule sends the ridge to the first lambda of its
# path at low events, a null model whose IBS is the same for every panel (?score_benchmark)
learner_of <- function(name) switch(name,
  l2cox = mlr3::lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min"),
  rsf   = mlr3::lrn("surv.ranger", num.trees = 500))

grid_selectors <- function(learner, cfg) {
  s <- default_selectors(prop = cfg$prop, learner = learner, q_grid = cfg$q_grid,
                         cv_args = list(inner_folds = cfg$inner_folds, B = cfg$B))
  for (q in cfg$q_fixed) s[[sprintf("survfcbf_q%d", q)]] <- local({
    qq <- q; bb <- cfg$B
    function(ti, st, X) survfcbf_select(ti, st, X, q = qq, B = bb)$selected
  })
  if (!is.null(cfg$selectors)) s <- s[intersect(cfg$selectors, names(s))]
  s
}

grid_cells <- function(cfg) {
  scn <- build_scenarios(cfg)
  cells <- expand.grid(s = seq_along(scn), learner = cfg$learners, rep = seq_len(cfg$R),
                       stringsAsFactors = FALSE)
  cells$scenario <- vapply(cells$s, function(i) scn[[i]]$name, character(1))
  cells
}

run_cell <- function(scn, learner_name, rep, cfg, master_seed) {
  if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
  seed <- as.integer(master_seed * 1000L + rep)
  sim <- simulate_dataset(utils::modifyList(scn$params, list(seed = seed)))
  task <- surv_task(sim$time, sim$status, sim$X)
  L <- learner_of(learner_name)
  tab <- benchmark_selectors(task, grid_selectors(L, cfg), L$clone(deep = TRUE),
                             outer_folds = cfg$outer, manifest = sim$manifest, seed = seed)
  data.frame(scenario = scn$name, map = scn$map, learner = learner_name, replicate = rep,
             events = sum(sim$status), tab, row.names = NULL, check.names = FALSE)
}

aggregate_raw <- function(raw) {
  metrics <- intersect(c("uno_c", "harrell_c", "ibs", "stability", "parsimony", "runtime",
                         "sensitivity", "fdr", "events"), names(raw))
  keys <- unique(raw[c("scenario", "map", "learner", "selector")])
  se <- function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    sub <- merge(k, raw)
    for (m in metrics) {
      x <- sub[[m]][is.finite(sub[[m]])]
      k[[paste0(m, "_mean")]] <- if (length(x)) mean(x) else NA_real_
      k[[paste0(m, "_se")]]   <- se(x)
    }
    k
  }))
}
