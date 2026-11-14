# TEST 1: does the instability MATTER? (the evidence gap). The benchmark's
# stability is already correlation-aware (stabm::stabilityIntersectionCount, theta=0.9), so
# survfcbf's low score means it picks genuinely UNRELATED genes across resamples -- not
# correlated substitutes. Open question: do those different genes nonetheless converge on
# the same BIOLOGY? Here each selection is mapped to the set of MSigDB pathways its genes
# touch, and stability is recomputed at the PATHWAY level and contrasted with the gene level.
#
# The control is what makes this test meaningful: pathways are large and overlapping, so
# even RANDOM gene sets touch common pathways by chance and would show spurious "functional
# stability". A random selector is therefore included, and only a functional stability
# clearly ABOVE the random control counts as convergence. SLURM array, coding-only.
suppressMessages({ library(survfcbf); library(farff); library(jsonlite); library(curl); library(stabm) })
bench <- Sys.getenv("BENCH_DIR", unset = ""); if (nzchar(bench)) setwd(bench)
source("realdata/load_bommert.R"); Sys.setenv(BOMMERT_DATA = file.path(bench, "realdata/Data"))
reg <- baseline_registry(); ent <- function(nm) Filter(function(e) e$name == nm, reg)[[1]]
cohorts <- names(BOMMERT_IDS)
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")); co <- cohorts[id]
outdir <- "results/realdata_functional"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(outdir, paste0(co, ".rds")); if (file.exists(out)) { cat(co, "skip\n"); quit(save = "no") }

d <- load_bommert(co); X <- d$X; nm <- colnames(X)
pc <- readRDS("realdata/prior_knowledge/protein_coding.rds")$ensembl
X <- X[, grepl("_rna$", nm) & sub("_rna$", "", nm) %in% pc, drop = FALSE]
ens <- sub("_rna$", "", colnames(X))
P   <- readRDS("realdata/prior_knowledge/ens2pathways.rds")
paths_all <- P$pathways; np <- length(paths_all)
gene_paths <- P$ens2path[ens]                          # per-column pathway membership
cat(sprintf("%s: n=%d ev=%d p=%d | %d genes with pathway annotation\n",
            co, length(d$time), sum(d$status), ncol(X), sum(lengths(gene_paths) > 0)))

K <- 10L; B <- 20L
## pathway "signature" of a gene selection: indices of pathways touched by >=1 selected gene
path_set <- function(idx) {
  p <- unique(unlist(gene_paths[idx], use.names = FALSE))
  match(p[!is.na(p)], paths_all)
}
selectors <- list(
  survfcbf_topk = function(t, s, Xs) survfcbf_select(t, s, Xs, q = 2L, top = K)$selected,
  variance      = function(t, s, Xs) cut_ranking(apply(Xs, 2L, stats::var), k = K),
  carss         = function(t, s, Xs) cut_ranking(abs(run_baseline(ent("carss"), t, s, Xs, prop = NULL, seed = 1L)$score), k = K),
  random        = function(t, s, Xs) sample.int(ncol(Xs), K)   # CONTROL: chance-level floor
)
G <- setNames(vector("list", length(selectors)), names(selectors))   # gene selections
Pa <- G                                                              # pathway signatures
set.seed(1)
for (b in seq_len(B)) {
  i <- sample(length(d$time), round(0.8 * length(d$time)))
  Xs <- X[i, , drop = FALSE]; ts <- d$time[i]; ss <- d$status[i]
  for (sl in names(selectors)) {
    g <- selectors[[sl]](ts, ss, Xs)
    G[[sl]][[b]]  <- g
    Pa[[sl]][[b]] <- path_set(g)
  }
}
res <- do.call(rbind, lapply(names(selectors), function(sl) data.frame(
  cohort = co, n = length(d$time), events = sum(d$status), p = ncol(X), selector = sl,
  stab_gene    = stabm::stabilityNogueira(G[[sl]],  p = ncol(X)),
  stab_pathway = stabm::stabilityNogueira(Pa[[sl]], p = np),
  paths_touched = mean(lengths(Pa[[sl]])),
  genes_annotated = mean(vapply(G[[sl]], function(g) sum(lengths(gene_paths[g]) > 0), numeric(1))),
  stringsAsFactors = FALSE)))
saveRDS(res, out)
print(res[, c("selector", "stab_gene", "stab_pathway", "paths_touched")], row.names = FALSE)
cat("FUNCTIONAL_DONE\n")
