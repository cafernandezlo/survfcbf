# experiments

A `targets` pipeline that exercises the installed `survfcbf` simulator. It depends on the package being installed
(`devtools::install()` from the repo root), not on the source tree.

`_targets.R` runs the default scenario, the censoring sweep (0.2 / 0.5 / 0.7) and
the low-variance-signal scenario through `run_replicates()`, assembles the
realised-diagnostics table (censoring fraction, concordance, realised
intra-cluster correlation, effective events) with Monte Carlo standard errors,
and sanity-checks the low-variance-signal scenario.

Run from this directory:

```r
targets::tar_make()
targets::tar_read(diagnostics_table)
targets::tar_read(low_variance_sanity)
```

Settings: `R = 50` replicates for this demonstration (the ficha production value
is `R = 200`). The `_targets/` store is git-ignored; the produced table is
written to `results/diagnostics_table.csv`.
