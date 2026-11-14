# Regenerate the results/ allowlist in experiments/benchmark/.gitignore.
#
# WHY THIS EXISTS. results/ used to be denied wholesale, so a fresh clone could not
# regenerate a single number of the Application Note -- which Oxford Bioinformatics treats as a
# submission blocker, not a convenience ("all data on which the conclusions given in the
# publication are based must be publicly available in stable public repositories").
# Tracking the whole 13 MB tree instead is the wrong trade: .rds files do not delta
# compress, so every pipeline re-run would append a full copy to history forever, and
# realdata_scores/ alone is 5 MB that feeds no published number.
#
# WHAT IT DOES. Runs paper/make_statistics.R, paper/make_appnote_figure.R and
# paper/make_supplementary_tables.R, the scripts that produce every number, the figure
# and the supplementary tables of the Application Note, with the read
# primitives instrumented, and keeps exactly the files under results/ that were actually
# opened. The list is therefore derived, not curated: if someone adds a statistic sourced
# from a new file, re-running this picks it up. Editing the .gitignore by hand instead is
# how the allowlist silently drifts out of sync with the scripts.
#
# The primitives are rebound INSIDE their namespaces. A globalenv override is not enough:
# make_statistics.R calls utils::read.csv() namespace-qualified, which bypasses it, and
# the resulting allowlist silently omits benchmark_realdata.csv and the ORA tables.
#
# Usage (from the repository root):
#   Rscript --no-init-file experiments/benchmark/derive_tracked_results.R
# Then verify the result matches, which is the check that actually matters:
#   git ls-files --others --exclude-standard experiments/benchmark/results/

ROOT <- normalizePath(Sys.getenv("SURVFCBF_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
RES  <- file.path(ROOT, "experiments/benchmark/results")
GI   <- file.path(ROOT, "experiments/benchmark/.gitignore")
stopifnot("run from the repository root (or set SURVFCBF_ROOT)" = dir.exists(RES))

READS <- new.env(parent = emptyenv()); READS$p <- character(0)
norm <- function(p) gsub("\\\\", "/", suppressWarnings(normalizePath(p, winslash = "/", mustWork = FALSE)))
note <- function(v) { if (is.character(v) && length(v)) READS$p <- c(READS$p, norm(v)); invisible() }
assign("note", note, envir = globalenv())

patched <- list()
patch <- function(ns, name, grab) {
  env <- asNamespace(ns); orig <- get(name, envir = env)
  was <- bindingIsLocked(name, env)
  if (was) unlockBinding(name, env)
  assign(name, grab(orig), envir = env)
  if (was) lockBinding(name, env)
  patched[[length(patched) + 1L]] <<- list(ns = ns, name = name, orig = orig, was = was)
}
unpatch_all <- function() for (p in rev(patched)) {
  env <- asNamespace(p$ns)
  if (p$was) unlockBinding(p$name, env)
  assign(p$name, p$orig, envir = env)
  if (p$was) lockBinding(p$name, env)
}

patch("utils", "read.table",  function(o) function(file, ...) { note(file); o(file, ...) })
patch("utils", "read.csv",    function(o) function(file, ...) { note(file); o(file, ...) })
patch("base",  "readRDS",     function(o) function(file, ...) { note(file); o(file, ...) })
patch("base",  "readLines",   function(o) function(con = stdin(), ...) { note(con); o(con, ...) })
patch("base",  "file.exists", function(o) function(...) { note(unlist(list(...))); o(...) })
patch("base",  "list.files",  function(o) function(path = ".", ...) {   # record what it RETURNS
  out <- o(path = path, ...); a <- list(...)
  if (length(out)) note(if (isTRUE(a$full.names)) out else file.path(path, out))
  out
})

run <- function(wd, script) {
  old <- getwd(); on.exit(setwd(old)); setwd(file.path(ROOT, wd))
  tryCatch({ source(script, local = new.env(parent = globalenv())); TRUE },
           error = function(e) { message("  ERROR in ", script, ": ", conditionMessage(e)); FALSE })
}
Sys.setenv(SURVFCBF_ROOT = ROOT)
ok1 <- run(".", "paper/make_statistics.R")
ok2 <- run(".", "paper/make_appnote_figure.R")
ok3 <- run(".", "paper/make_supplementary_tables.R")
unpatch_all()
if (!ok1 || !ok2 || !ok3) stop("a consuming script failed; the allowlist would be incomplete")

hit <- unique(READS$p)
hit <- hit[startsWith(hit, norm(RES))]
hit <- hit[file.exists(hit) & !dir.exists(hit)]
keep <- sort(unique(sub(paste0("^", norm(RES), "/"), "", hit)))

hdr <- c(
  "# Benchmark outputs. Deny by default, then re-include EXACTLY the files that",
  "# paper/make_statistics.R, paper/make_appnote_figure.R and",
  "# paper/make_supplementary_tables.R read, so a fresh clone can regenerate every number,",
  "# the figure and the supplementary tables of the Application Note with no cluster and",
  "# no downloads.",
  "#",
  "# Derived, not curated: regenerate with",
  "#   Rscript --no-init-file experiments/benchmark/derive_tracked_results.R",
  "# rather than editing by hand. Note git cannot re-include a file whose parent",
  "# directory is excluded, which is why each directory is un-excluded before its",
  "# contents are filtered.",
  "#",
  "# Deliberately NOT tracked: realdata_scores/ (5.0 MB, feeds only realdata_gsea/gsea.rds,",
  "# which is tracked), realdata_signatures*/backgrounds.rds, gsea_full.csv and the other",
  "# intermediates. The full tree goes to the Zenodo deposit instead.",
  "results/*")
body <- sprintf("!results/%s", keep[!grepl("/", keep)])
for (d in sort(unique(dirname(keep[grepl("/", keep)])))) {
  fs <- sort(basename(keep[dirname(keep) == d]))
  body <- c(body, "", sprintf("!results/%s/", d), sprintf("results/%s/*", d),
            sprintf("!results/%s/%s", d, fs))
}
## Rewrite our own block and keep every unrelated rule. The filter must drop EVERY
## results/ line, not just the blanket deny: on a second run, leftover `results/<dir>/*`
## lines would otherwise be appended AFTER the negations and silently re-exclude
## everything the allowlist just re-included.
old  <- readLines(GI, warn = FALSE)
t    <- trimws(old)
mine <- startsWith(t, "#") | startsWith(t, "results/") | startsWith(t, "!results/") | !nzchar(t)
rest <- old[!mine]
writeLines(c(hdr, body, "", rest), GI)

cat(sprintf("allowlist: %d files, %.2f MB on disk\n", length(keep),
            sum(file.size(file.path(RES, keep))) / 1024^2))
cat("check with: git ls-files --others --exclude-standard experiments/benchmark/results/\n")
