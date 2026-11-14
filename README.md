# survfcbf

Fast correlation-based filter (FCBF) feature selection for right-censored survival data, with a
leakage-safe harness to evaluate it and any other survival filter.

FCBF (Yu and Liu, 2003) keeps a feature only if no more relevant feature already explains it, using
an approximate Markov blanket. `survfcbf` applies that rule to censored outcomes through the
martingale residuals of a null Cox model. The package also ships:

- a nested cross-validation harness in which everything learned from data stays inside the training
  fold,
- a survival simulator with known ground truth, scored at the level of correlated clusters,
- uniform wrappers for the standard baselines: the variance filter, correlation-adjusted regression
  scores, mutual-information rankers, penalised Cox and random survival forests.

**When to use it.** `survfcbf` is for the case where a small, non-redundant, cheaply obtained panel is
what you need. It is not the best choice for discrimination whatever the panel size, nor for a
ranking that stays stable under resampling. On real expression data the rankers the package wraps do
better on both counts, and the accompanying Application Note reports these comparisons in full.

## Installation

```r
# install.packages("remotes")
remotes::install_github("cafernandezlo/survfcbf", build_vignettes = TRUE)
```

The filter and the simulator need only `survival`. The evaluation harness and the baselines use
optional packages. `mlr3proba` and `mlr3extralearners` are not on CRAN, so install them from the
mlr-org r-universe, keeping CRAN as a second repository for their dependencies:

```r
install.packages(c("mlr3", "glmnet", "praznik", "stabm", "carSurv", "ranger", "grpreg", "lgr"))
install.packages(c("mlr3proba", "mlr3extralearners"),
                 repos = c("https://mlr-org.r-universe.dev", "https://cloud.r-project.org"))
```

The package tests compare the filter with selections of the reference implementation, the `FCBF`
package, frozen in the test file, so they do not need it. Only the full validation,
`experiments/benchmark/realdata/fcbf_crossvalidation.R`, runs `FCBF` itself. It left Bioconductor
with release 3.20, so install it from GitHub:

```r
remotes::install_github("lubianat/FCBF")
```

## Quick start

Select a panel of ten features from a simulated dataset:

```r
library(survfcbf)

sim <- simulate_dataset(default_params(
  p = 300L, n = 400L, k_drivers = 3L, cluster_size = 5L,
  n_distractor_subblocks = 2L, distractor_size = 8L,
  censor_fraction = 0.5, target_C = 0.75,
  seed = 1L, calib_n_presample = 800L
))

sel <- survfcbf_select(sim$time, sim$status, sim$X, q = 2L, top = 10L)
sel$selected                                 # column indices of the panel
head(sel$trace$removed)                      # what was pruned, and by which feature
score_selection_clusterwise(sel$selected, sim$manifest)$sensitivity
```

`top` fixes the panel size and is the recommended mode: it is deterministic and needs no permutations.
Without it the filter sets its own size from a permutation null and can return empty panels when
events are scarce. `q` is the number of bins. The results reported for the package use `q = 2` with
`top`, while the default `q = NULL` follows Bommert's rule, which gives 10 from 30 samples up.

Evaluate selectors under nested cross-validation:

```r
library(mlr3); library(mlr3proba); library(mlr3extralearners)

task <- surv_task(sim$time, sim$status, sim$X)
selectors <- default_selectors(prop = 0.1, topk = 10L)[c("survfcbf_topk", "variance", "mim_surv")]
learner <- lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")

benchmark_selectors(task, selectors, learner, outer_folds = 3L, seed = 1L)
```

Keep `s = "lambda.min"`. At its default `lambda.1se`, with far more features than events, the ridge Cox
learner collapses to the null model and every panel gets the same integrated Brier score. The
`n_null_fit` column counts the folds where this happened.

For a guided tour, including how to read the output, see the vignette:

```r
vignette("survfcbf", package = "survfcbf")
```

## Main functions

| Task | Functions |
|---|---|
| Select features | `survfcbf_select()`, `survfcbf_cv_select()`, `fcbf_select()` |
| Evaluate selectors | `benchmark_selectors()`, `nested_cv_select()`, `score_benchmark()`, `benchmark_stability()`, `bayes_signed_rank()` |
| Baselines | `default_selectors()`, `baseline_registry()`, `run_baseline()`, `variance_select()`, `carss_select()`, `mim_surv_select()`, `mrmr_surv_select()`, `cox_glmnet_select()`, `rsf_select()` |
| Simulate | `simulate_dataset()`, `default_params()`, `run_replicates()`, `score_selection_clusterwise()` |
| Building blocks | `symmetrical_uncertainty()`, `discretize_quantile()`, `survival_target()`, `delta_permutation()` |

## Reproducing the benchmark

The pipeline lives in `experiments/benchmark` and is not part of the installed package. The result
tables are tracked under `experiments/benchmark/results`.

```bash
cd experiments/benchmark
Rscript --no-init-file realdata/run_all_local.R      # analyses that run locally
bash realdata/submit_all_cesga.sh --with-benchmark   # SLURM arrays, about 12-14 h per cohort
```

The SLURM half is written for our cluster. Its job scripts look for this directory and the R
package library under `SURVFCBF_CLUSTER_ROOT` (default `$STORE2/survfcbf`), as
`experiments/benchmark` and `Rlib`, and `SURVFCBF_RLIB` overrides the library alone. The container
and module names are site-specific. The data are not in the repository:

- **Eleven TCGA-derived cohorts**, prepared by Herrmann et al. (2020) and selected by Bommert et al.
  (2022): `realdata/load_bommert.R` downloads them from OpenML and caches them in `BOMMERT_DATA`.
- **METABRIC**, the complete cBioPortal profile `brca_metabric_mrna`: `realdata/load_metabric.R` reads
  the expression snapshot and the clinical table from `METABRIC_DIR` and stops if either SHA-256 does
  not match. A later download may not match byte for byte. The DOI of the archived snapshot will be
  added here at release.

The figures quoted in the Application Note come from a single file and are checked against the text:

```bash
Rscript --no-init-file paper/make_statistics.R
Rscript --no-init-file paper/make_appnote_figure.R
Rscript --no-init-file paper/check_manuscript_numbers.R
```

## License

The package is GPL-3. METABRIC data are under the ODC Open Database License as distributed by
cBioPortal.

## Citation

Application Note in preparation. Until then, please cite the repository.
