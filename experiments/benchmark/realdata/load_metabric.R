# Loader for the METABRIC external cohort (external validation +
# the real-data test of the events-dependency mechanism). cBioPortal-derived local .rds:
# expression (gene symbols x MB samples, log2 Illumina microarray), clinical (OS). Returns
# list(time, status, X) with X = samples x genes, OS endpoint, NA/zero-variance genes
# dropped. Path via the METABRIC_DIR environment variable, which must point at the folder
# holding both snapshots. assert_scale guards against a wrong (non-log2) load. No cross-platform transfer is
# done -- METABRIC is benchmarked on its own, so the microarray-vs-RNAseq batch is a
# non-issue.
#
# Expression is the COMPLETE profile brca_metabric_mrna (Illumina HT-12 v3, 19,850 genes x
# 1,980 samples, 19,839 genes after the NA/variance filter), the snapshot
# metabric_expression_full.rds downloaded from the cBioPortal REST API (portal v7.2.0) on
# 2026-10-25. Its provenance is described in the README of the METABRIC Zenodo record,
# doi:10.5281/zenodo.24473297, which holds both snapshots.
# Both the expression snapshot and the clinical table are SHA-256 checked, so a swapped
# input stops the run instead of changing the numbers.
load_metabric <- function(dir = Sys.getenv("METABRIC_DIR", ""),
                          file = "metabric_expression_full.rds",
                          sha256 = "4acb9a0027af72fb2a85c1cd1e8eee1e2d82f21367f647c6fe4ff8078480ee8c",
                          clinical = "metabric_clinical.rds",
                          clinical_sha256 = "5572d56677e8d59bad60fc0055cd0e06a221aeb22cd0ef09a1c016070eee4bb7") {
  if (!nzchar(dir))
    stop("set METABRIC_DIR to the folder holding ", file, " and ", clinical,
         " (both are in the Zenodo record doi:10.5281/zenodo.24473297)", call. = FALSE)
  if (!dir.exists(dir)) stop("METABRIC_DIR does not exist: ", dir, call. = FALSE)
  checked <- function(name, expected) {
    f <- file.path(dir, name)
    if (!file.exists(f)) stop("METABRIC input not found: ", f, call. = FALSE)
    if (!identical(unname(tools::sha256sum(f)), expected))
      stop("METABRIC input SHA-256 mismatch (expected ", expected, "): ", f, call. = FALSE)
    readRDS(f)
  }
  e  <- checked(file, sha256)                                 # genes(symbol) x samples
  cl <- checked(clinical, clinical_sha256)
  common <- intersect(colnames(e), cl$patientId)
  sub <- cl[match(common, cl$patientId), ]
  time   <- suppressWarnings(as.numeric(sub$OS_MONTHS))
  status <- as.integer(grepl("DECEASED|^1", sub$OS_STATUS))   # OS endpoint (1 = death)
  X <- t(as.matrix(e[, common, drop = FALSE]))                # samples x genes
  keep_gene <- colSums(is.na(X)) == 0L & apply(X, 2L, stats::var) > 0
  X <- X[, keep_gene, drop = FALSE]
  ok <- is.finite(time) & time > 0 & !is.na(status)
  X <- X[ok, , drop = FALSE]; time <- time[ok]; status <- status[ok]
  stopifnot("scale mismatch: expected log2 microarray (max < 25)" = max(X, na.rm = TRUE) < 25)
  list(time = time, status = status, X = X)
}
