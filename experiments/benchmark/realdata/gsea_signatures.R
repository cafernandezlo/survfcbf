# GSEA of each selector's full per-gene relevance ranking (the
# ranking-based interpretability test that complements the top-K ORA). For each cohort x
# selector, rank the coding genes by the selector's relevance score (Ensembl -> symbol)
# and run fgsea against MSigDB Hallmarks + oncogenic C6. Relevance scores are non-negative
# magnitudes (SU, variance, MI; |CAR| for carss), so a one-sided GSEA (scoreType="pos")
# asks the mechanistic question: are pathway genes concentrated among the highest-relevance
# genes? fgsea's padj is BH within each cohort x selector ranking. Run AFTER extract_scores.R
# and pulling results/realdata_scores/.
bench <- Sys.getenv("BENCH_DIR", "")
if (nzchar(bench) && dir.exists(bench)) setwd(bench)
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi); library(fgsea) })
set.seed(1)

SCOREDIR <- Sys.getenv("RD_SCORES_DIR", "results/realdata_scores")
PK  <- "realdata/prior_knowledge"
OUT <- "results/realdata_gsea"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msig <- readRDS(file.path(PK, "msigdb_H_C6.rds"))
sets <- split(msig$sets$gene_symbol, msig$sets$gs_name)
selectors <- c("survfcbf", "variance", "carss", "mim_surv")

files <- list.files(SCOREDIR, pattern = "rds$", full.names = TRUE)
res <- list(); maps <- list()
## SU at q = 2 takes few distinct values, so the survfcbf ranking is full of ties. They
## are broken as survfcbf_select() breaks them, by the absolute Spearman correlation with
## the martingale residuals, added in a step smaller than half the smallest gap between
## distinct SU values so the SU order itself never changes.
tie_broken <- function(su, tb) {
  u <- sort(unique(su)); gap <- if (length(u) > 1L) min(diff(u)) else 1
  su + tb * gap / 2
}
for (f in files) {
  co <- sub("\\.rds$", "", basename(f)); sc <- readRDS(f)
  distinct_su <- length(unique(round(sc$survfcbf, 12)))           # before the tie-break
  if (!is.null(sc$survfcbf_tiebreak)) sc$survfcbf <- tie_broken(sc$survfcbf, sc$survfcbf_tiebreak)
  ens <- sub("_rna$", "", sc$feature)
  sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  ok <- !is.na(sym) & !duplicated(sym)
  maps[[co]] <- data.frame(cohort = co, coding_features = length(ens), no_symbol = sum(is.na(sym)),
                           duplicated_symbol = sum(!is.na(sym) & duplicated(sym)), ranked = sum(ok),
                           distinct_su = distinct_su)
  for (sl in selectors) {
    v <- abs(sc[[sl]][ok]); names(v) <- sym[ok]; v <- v[is.finite(v)]
    fg <- as.data.frame(suppressWarnings(fgsea(pathways = sets, stats = sort(v, decreasing = TRUE),
                                               minSize = 5, maxSize = 500, scoreType = "pos")))
    if (nrow(fg))
      res[[length(res) + 1L]] <- data.frame(cohort = co, selector = sl, pathway = fg$pathway,
        NES = fg$NES, pval = fg$pval, padj = fg$padj, size = fg$size, stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, res)
utils::write.csv(R, file.path(OUT, "gsea_full.csv"), row.names = FALSE)
saveRDS(R, file.path(OUT, "gsea.rds"))
utils::write.csv(do.call(rbind, maps), file.path(OUT, "mapping.csv"), row.names = FALSE)

cat("=== GSEA: pathways enriched (padj<0.05) per selector, K-free ranking test ===\n")
summ <- do.call(rbind, lapply(selectors, function(sl) {
  rr <- R[R$selector == sl, ]
  hits <- tapply(rr$padj < 0.05, rr$cohort, function(x) sum(x, na.rm = TRUE))
  data.frame(selector = sl, cohorts_enriched = sum(hits >= 1), mean_paths = round(mean(hits), 1),
             total_sig = sum(rr$padj < 0.05, na.rm = TRUE), stringsAsFactors = FALSE)
}))
print(summ, row.names = FALSE)
cat("\n=== survfcbf top enriched pathways (padj<0.05) ===\n")
tp <- R[R$selector == "survfcbf" & !is.na(R$padj) & R$padj < 0.05, c("cohort", "pathway", "NES", "padj")]
if (nrow(tp)) print(head(tp[order(tp$padj), ], 12), row.names = FALSE) else cat("  none\n")
cat(sprintf("\ntables -> %s\nGSEA_DONE\n", OUT))
