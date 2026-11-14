# ORA of the module-FCBF selected modules (module pivot, interpretability
# pillar). Each selected co-expression module is a gene set; a module "enriches" if it hits
# >=1 MSigDB Hallmark/C6 pathway at BH q<0.05 (per-module BH over the pathways). Reports,
# per relevance variant (SU-ranked vs CAR-ranked modules), how many cohorts have >=1
# enriched module and the mean # enriched modules -- the contrast with the gene-level
# top-10 that enriched 0/11. Local; consumes results/realdata_modules_v3/*.rds + the cached
# MSigDB. rna genes are Ensembl -> symbol via org.Hs.eg.db; background = the cohort's genes.
bench <- Sys.getenv("BENCH_DIR", ""); if (nzchar(bench) && dir.exists(bench)) setwd(bench)
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
msig <- readRDS("realdata/prior_knowledge/msigdb_H_C6.rds")
sets <- split(msig$sets$gene_symbol, msig$sets$gs_name)
files <- list.files("results/realdata_modules_v3", pattern = "rds$", full.names = TRUE)
Q <- 0.05; MIN_SET <- 5L

map_syms <- function(ens) { s <- mapIds(org.Hs.eg.db, keys = unique(ens), keytype = "ENSEMBL",
                                        column = "SYMBOL", multiVals = "first"); unique(s[!is.na(s)]) }
## does a module (symbol set) hit >=1 pathway at BH q<0.05, given a background?
module_enriches <- function(S, B) {
  N <- length(B); nS <- length(S); if (nS < 2L) return(FALSE)
  p <- vapply(sets, function(P) {
    Pb <- intersect(P, B); m <- length(Pb); if (m < MIN_SET) return(NA_real_)
    stats::phyper(length(intersect(S, Pb)) - 1L, m, N - m, nS, lower.tail = FALSE)
  }, numeric(1))
  q <- stats::p.adjust(p[!is.na(p)], "BH"); any(q < Q)
}

rows <- list()
for (f in files) {
  r <- readRDS(f); co <- r$summary$cohort
  B <- map_syms(r$background)
  for (variant in c("mods_su", "mods_car")) {
    enr <- vapply(r[[variant]], function(g) module_enriches(map_syms(g), B), logical(1))
    rows[[length(rows) + 1L]] <- data.frame(cohort = co, variant = sub("mods_", "", variant),
      n_modules_sel = length(r[[variant]]), n_enriched = sum(enr),
      any_enriched = any(enr), stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, rows)
OUT <- "results/realdata_modules_v3"; utils::write.csv(R, file.path(OUT, "module_enrichment.csv"), row.names = FALSE)
cat("=== module enrichment (H+C6, per-module BH q<0.05) ===\n")
for (v in c("su", "car")) {
  rr <- R[R$variant == v, ]
  cat(sprintf("  module-%-3s: %d/%d cohorts with >=1 enriched module | mean %.1f enriched of 10\n",
              toupper(v), sum(rr$any_enriched), nrow(rr), mean(rr$n_enriched)))
}
cat("(gene-level top-10 ORA was 0/11 cohorts)\n")
print(R[R$variant == "su", c("cohort", "n_enriched", "any_enriched")], row.names = FALSE)
cat("MODULE_ENRICH_DONE\n")
