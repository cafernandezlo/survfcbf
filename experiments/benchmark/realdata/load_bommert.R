# Loader for the 11 high-dimensional gene-expression survival datasets of
# Bommert et al. (2022, Briefings in Bioinformatics). Each dataset
# is a set of OpenML ARFF "parts" (public, TCGA-derived, de-identified); we download
# and cache them, then return (time, status, X) where X is the rna + mirna expression
# features -- exactly Bommert's set: their grep("rna") catches both, e.g. PAAD has
# 22348 rna + 612 mirna = 22960 = their reported p. Other modalities (cnv, mutation,
# clinical, barcode) are dropped. Data is never committed; it lives in the cache dir.
#
# Deps: jsonlite, curl, farff. Cache dir from env BOMMERT_DATA (else "Data" in cwd).
# Usage:  source("load_bommert.R"); d <- load_bommert("PAAD")
#         surv_task(d$time, d$status, d$X)   # feeds the same harness

BOMMERT_IDS <- list(                                  # OpenML data IDs per cohort
  BLCA = c(42279, 42280), BRCA = c(42323, 42324, 42325), HNSC = c(42285, 42286),
  KIRC = c(42287, 42288), LGG = c(42293, 42294), LUAD = c(42297, 42298),
  LUSC = c(42299, 42300), OV = c(42301, 42302), PAAD = c(42303, 42304),
  SKCM = c(42311, 42312), STAD = c(42313, 42314)
)

.bommert_cache <- function(cache_dir) {
  if (is.null(cache_dir)) cache_dir <- Sys.getenv("BOMMERT_DATA", "Data")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

.arff_url <- function(id)
  jsonlite::fromJSON(sprintf("https://www.openml.org/api/v1/json/data/%d", id))$data_set_description$url

.robust_dl <- function(url, dest, tries = 6) {            # OpenML drops large transfers; retry, atomic rename
  for (t in seq_len(tries)) {
    ok <- tryCatch({
      curl::curl_download(url, paste0(dest, ".part"), quiet = TRUE, mode = "wb")
      file.rename(paste0(dest, ".part"), dest)
      file.size(dest) > 1e6
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(invisible(TRUE))
    Sys.sleep(3)
  }
  stop("download failed after ", tries, " tries: ", url, call. = FALSE)
}

## MD5 of every ARFF part as cached for the published runs (CESGA, 2026-11-19). A part whose
## checksum differs stops the loader, so a silently re-uploaded or truncated OpenML file
## cannot change a result. BOMMERT_MD5 points elsewhere, or to "" to skip the check.
.bommert_md5 <- function(f) {
  manifest <- Sys.getenv("BOMMERT_MD5", "realdata/openml_md5.csv")
  if (!nzchar(manifest) || !file.exists(manifest)) return(invisible(NULL))
  m <- utils::read.csv(manifest, stringsAsFactors = FALSE)
  want <- m$md5[m$file == basename(f)]
  if (length(want) != 1L) stop("no checksum for ", basename(f), " in ", manifest, call. = FALSE)
  got <- unname(tools::md5sum(f))
  if (!identical(got, want)) stop("checksum mismatch for ", basename(f), ": ", got, " != ", want, call. = FALSE)
  invisible(NULL)
}

#' Download (and cache) the OpenML ARFF parts for one Bommert cohort; returns the file paths.
download_bommert <- function(name, cache_dir = NULL) {
  cache_dir <- .bommert_cache(cache_dir)
  ids <- BOMMERT_IDS[[name]]
  if (is.null(ids)) stop("unknown dataset '", name, "'; one of: ", paste(names(BOMMERT_IDS), collapse = ", "))
  vapply(ids, function(id) {                              # name by OpenML id, not URL basename:
    f <- file.path(cache_dir, sprintf("%s_%d.arff", name, id))  # OpenML reuses basenames across ids
    if (!file.exists(f)) .robust_dl(.arff_url(id), f)     # (BRCA 42324 and 42325 are both "BRCA-Part2.arff",
    .bommert_md5(f)                                       # which silently overwrote part 42324) -- and lets a
    f                                                     # cached file be found without re-hitting the API
  }, character(1))

}

#' Load one Bommert cohort into (time, status, X) with X = rna + mirna expression.
load_bommert <- function(name, cache_dir = NULL) {
  parts <- lapply(download_bommert(name, cache_dir), function(f) suppressMessages(farff::readARFF(f)))
  dat <- do.call(cbind, parts)
  keep <- grep("rna$", colnames(dat), value = TRUE)     # rna + mirna (both end in "rna"); drops cnv/mutation/clinical
  X <- as.matrix(dat[, keep]); storage.mode(X) <- "double"
  list(time = as.numeric(dat[["time"]]), status = as.integer(dat[["status"]]), X = X, name = name,
       barcode = as.character(dat[["bcr_patient_barcode"]]))
}
