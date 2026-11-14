# Single entry point for the LOCAL half of the real-data stage. The
# compute-heavy halves are SLURM arrays on CESGA (see submit_all_cesga.sh); everything that
# analyses their outputs, plus the two self-contained experiments, runs here in dependency
# order. Each stage declares the inputs it needs, so a missing CESGA download is reported
# as a skip with the reason instead of a mid-run error, and the pipeline can be re-run
# incrementally as results arrive.
#
# Usage (from experiments/benchmark):
#   Rscript --no-init-file realdata/run_all_local.R           # every ready stage
#   RD_ONLY=cross_dataset,gsea Rscript --no-init-file realdata/run_all_local.R
#   RD_FORCE=1 Rscript ...                                    # re-run even if outputs exist
#
# sessionInfo() is written to results/session_info_local.txt on every run.
setwd(Sys.getenv("BENCH_DIR", unset = "."))
t_start <- proc.time()[["elapsed"]]
FORCE <- nzchar(Sys.getenv("RD_FORCE", ""))
ONLY  <- setdiff(strsplit(Sys.getenv("RD_ONLY", ""), ",")[[1]], "")

## stage: name, the script to run, inputs that must exist, and a representative output
stages <- list(
  list(name = "prior_knowledge", script = "realdata/cache_prior_knowledge.R",
       needs = character(), out = "realdata/prior_knowledge/msigdb_H_C6.rds",
       note = "one-time, needs internet (MSigDB via msigdbr, intOGen, HGNC)"),
  list(name = "aggregate", script = "realdata/aggregate_realdata.R",
       needs = "results/realdata/cells", out = "results/benchmark_realdata.csv",
       note = "needs the 11 per-cohort benchmark CSVs from CESGA"),
  list(name = "cross_dataset", script = "realdata/analyze_cross_dataset.R",
       needs = "results/benchmark_realdata.csv", out = "results/realdata_crossdataset/crossdataset.rds",
       note = "Friedman/Nemenyi CD + Bayesian signed-rank"),
  list(name = "enrichment", script = "realdata/enrich_signatures.R",
       needs = c("results/realdata_signatures/signatures.csv", "realdata/prior_knowledge/msigdb_H_C6.rds"),
       out = "results/realdata_enrichment/enrichment.rds",
       note = "ORA of the top-K signatures (set SIGDIR/ENRICH_OUT for the coding-only variant)"),
  list(name = "gsea", script = "realdata/gsea_signatures.R",
       needs = "results/realdata_scores", out = "results/realdata_gsea/gsea.rds",
       note = "GSEA of the full relevance rankings"),
  list(name = "module_enrichment", script = "realdata/module_enrich.R",
       needs = "results/realdata_modules_v3", out = "results/realdata_modules_v3/module_enrichment.csv",
       note = "ORA of the selected co-expression modules"),
  list(name = "cardinality_inference", script = "realdata/cardinality_inference.R",
       needs = "results/realdata_cardinality", out = "results/realdata_cardinality/cardinality_inference.csv",
       note = "Bayesian signed-rank on the cardinality curve"),
  list(name = "runtime", script = "realdata/runtime_scaling.R",
       needs = character(), out = "results/realdata_runtime/runtime_scaling.csv",
       note = "self-contained (simulated data)"),
  list(name = "fcbf_crossvalidation", script = "realdata/fcbf_crossvalidation.R",
       needs = character(), out = "results/validation/fcbf_crossvalidation.csv",
       note = "self-contained (needs FCBF from GitHub, see the script header)"),
  list(name = "metabric_curve", script = "realdata/metabric_event_curve.R",
       needs = character(), out = "results/realdata_metabric/event_curve.csv",
       note = "self-contained (local METABRIC via METABRIC_DIR)")
)

run_stage <- function(s) {
  if (length(ONLY) && !(s$name %in% ONLY)) return("filtered")
  missing <- s$needs[!file.exists(s$needs)]
  if (length(missing)) {
    cat(sprintf("  SKIP  %-22s missing: %s\n", s$name, paste(basename(missing), collapse = ", ")))
    return("skipped")
  }
  if (!FORCE && file.exists(s$out)) {
    cat(sprintf("  DONE  %-22s (output exists; RD_FORCE=1 to redo)\n", s$name))
    return("cached")
  }
  cat(sprintf("  RUN   %-22s %s\n", s$name, s$script))
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({ callr_ok <- system2(file.path(R.home("bin"), "Rscript"),
                     c("--no-init-file", shQuote(s$script)), stdout = "", stderr = "")
                   callr_ok == 0L },
                 error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("  %s  %-22s %.1f min\n", if (ok) "OK   " else "FAIL ", s$name,
              (proc.time()[["elapsed"]] - t0) / 60))
  if (ok) "ran" else "failed"
}

cat("=== survfcbf real-data pipeline (local stages) ===\n")
cat(sprintf("working dir: %s\n\n", getwd()))
status <- vapply(stages, run_stage, character(1))

dir.create("results", showWarnings = FALSE)
## Each stage runs as a subprocess, so the parent sessionInfo() alone records almost
## nothing. Pin the actual analysis stack by version explicitly alongside it.
stack <- c("survfcbf", "mlr3", "mlr3proba", "mlr3extralearners", "stabm", "glmnet",
           "carSurv", "praznik", "ranger", "grpreg", "survival", "farff", "jsonlite",
           "curl", "fgsea", "msigdbr", "AnnotationDbi", "org.Hs.eg.db")
vers <- vapply(stack, function(p)
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"),
  character(1))
writeLines(c(paste("generated:", format(Sys.time())), "",
             "== analysis stack (stages run as subprocesses; versions pinned here) ==",
             sprintf("  %-20s %s", stack, vers), "",
             "== parent session ==",
             capture.output(utils::sessionInfo())), "results/session_info_local.txt")

cat(sprintf("\n=== summary (%.1f min) ===\n", (proc.time()[["elapsed"]] - t_start) / 60))
print(data.frame(stage = vapply(stages, `[[`, character(1), "name"), status = status), row.names = FALSE)
skipped <- vapply(stages, `[[`, character(1), "name")[status == "skipped"]
if (length(skipped)) {
  cat("\nStages skipped for missing inputs (fetch them from CESGA, see submit_all_cesga.sh):\n")
  for (s in stages[status == "skipped"]) cat(sprintf("  %-22s %s\n", s$name, s$note))
}
cat("\nsessionInfo -> results/session_info_local.txt\n")
if (any(status == "failed")) quit(status = 1L)
