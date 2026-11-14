# Simulation benchmark grid

A `targets` pipeline that runs the selector sweep over the locked grid design and
produces the results table: per scenario x learner x selector, the discrimination
(Uno's C, Harrell's C, integrated Brier), stability (SMA-Count), parsimony,
selection runtime and cluster-level recovery (sensitivity / FDR), each with a Monte
Carlo standard error, tagged by which figure (`map`) each scenario feeds.

Run from this directory:

```r
targets::tar_make()
```

It depends on the installed `survfcbf` package and the mlr3 stack
(`mlr3`, `mlr3proba`, `mlr3extralearners`, `stabm`). After changing the package,
reinstall it (`devtools::install(quick = TRUE)`) so the pipeline picks it up.

## Tiers

`.tier` at the top of `_targets.R` selects the configuration:

- `"demo"` (default): a tiny slice (2 event levels, L2-Cox, R = 1, a fast selector
  subset) that validates the pipeline locally in a few minutes.
- `"production"`: the full grid for HPC -- reference clustered-with-
  redundancy, axes effective-events (n) x censoring x dimensionality p, OAT sweeps
  (redundancy off, PH vs AFT), all twelve `default_selectors` (survfcbf, the two
  inner-CV-tuned variants, the nine baselines) plus fixed-q survfcbf for the
  optimal-q figure, learners L2-Cox + RSF, R = 200.

## Maps (the headline figures)

- `A_events` / `A_censor`: recovery + prediction vs effective events. The operating
  threshold (where survfcbf's FDR / parsimony advantage emerges) and the safe
  failure below it.
- `B_dim`: cardinality + FDR vs `p`. survfcbf flat (scale-invariant) where the
  ranking baselines' FDR grows with `p`.
- the fixed-q survfcbf selectors give the optimal-q-vs-events figure (map C).
- `OAT`: redundancy off (where predominance has nothing to prune) and the AFT
  non-PH target (martingale-transform robustness).

## Notes

- **Compute.** The tuned survfcbf variants run a nested inner CV, so the production
  grid is firmly HPC-scale (the demo's 2-scenario slice already takes minutes).
  Run on the cluster with a `batchtools` backend; budget resources per cell.
- **Operating regime.** survfcbf's recovery tracks effective events per fold, not
  `n` or `p`. At realistic effect sizes (target C = 0.70) the events threshold is
  higher, so the low-`n` cells sit below it and survfcbf under-recovers (its
  self-cardinality config can pick too few features for a Cox model, recorded as a
  failed fold); the high-`n` cells show the crossover. Read the recovery columns
  alongside the discrimination metrics.
