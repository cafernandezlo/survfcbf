# Biological validation of the full-cohort signatures (the enrichment
# pillar). Three read-outs, cardinality-matched (K=10 primary; 5/20 sensitivity),
# comparing survfcbf_topk against the top-K of the ranking baselines:
#   (1) ORA over MSigDB Hallmarks + oncogenic C6 (hypergeometric, correct measured-gene
#       background per cohort), BH-corrected globally within the family;
#   (2) overlap with intOGen cancer drivers, cohort-specific (TCGA type -> intOGen type)
#       and global, hypergeometric;
#   (3) pathway diversity (distinct sets a signature touches).
# All local: org.Hs.eg.db for the Ensembl->symbol map, cached MSigDB/intOGen snapshots.
# No external service is contacted -> privacy-safe and reproducible. rna features
# are Ensembl (ENSG..._rna); mirna (..._mirna) are not in the coding-gene sets, so they
# are excluded from ORA/overlap and counted separately. Run AFTER extract_signatures.R
# and pulling results/realdata_signatures/.

bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench) && dir.exists(bench)) setwd(bench)
suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
set.seed(1)

SIGDIR <- Sys.getenv("SIGDIR", "results/realdata_signatures")       # coding variant overrides via env
PK     <- "realdata/prior_knowledge"
OUT    <- Sys.getenv("ENRICH_OUT", "results/realdata_enrichment"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
K_FOCUS <- 10L; Q_SIG <- 0.05; MIN_SET <- 5L

sig  <- utils::read.csv(file.path(SIGDIR, "signatures.csv"), stringsAsFactors = FALSE)
bg   <- readRDS(file.path(SIGDIR, "backgrounds.rds"))
msig <- readRDS(file.path(PK, "msigdb_H_C6.rds"))
intg <- readRDS(file.path(PK, "intogen_drivers.rds"))

## ---- ID mapping: rna=Ensembl->symbol; mirna excluded from ORA -------------------
is_rna <- function(f) grepl("_rna$",  f)
is_mir <- function(f) grepl("_mirna$", f)
ens_of <- function(f) sub("_rna$", "", f[is_rna(f)])
all_ens <- unique(unlist(lapply(bg, ens_of)))
ens2sym <- mapIds(org.Hs.eg.db, keys = all_ens, keytype = "ENSEMBL",
                  column = "SYMBOL", multiVals = "first")
to_syms <- function(features) { s <- ens2sym[ens_of(features)]; unique(s[!is.na(s)]) }

bg_sym <- lapply(bg, to_syms)                                   # measured universe per cohort
loss <- t(sapply(names(bg), function(co) c(
  rna_features = sum(is_rna(bg[[co]])), mapped_symbols = length(bg_sym[[co]]))))
cat("=== ID mapping (rna Ensembl -> symbol), losses per cohort ===\n"); print(loss)

sets <- split(msig$sets$gene_symbol, msig$sets$gs_name)         # MSigDB symbol lists
first <- !duplicated(msig$sets$gs_name)
set_coll <- setNames(msig$sets$collection[first], msig$sets$gs_name[first])

## ---- per-signature ORA + driver overlap ---------------------------------------
groups <- unique(sig[, c("cohort", "selector", "K")])
ora_rows <- list(); drv_rows <- list()
for (i in seq_len(nrow(groups))) {
  co <- groups$cohort[i]; sel <- groups$selector[i]; K <- groups$K[i]
  feats <- sig$feature[sig$cohort == co & sig$selector == sel & sig$K == K]
  S <- to_syms(feats); nmir <- sum(is_mir(feats))
  B <- bg_sym[[co]]; N <- length(B); nS <- length(S)
  for (gs in names(sets)) {
    Pb <- intersect(sets[[gs]], B); m <- length(Pb)
    if (m < MIN_SET) next
    k <- length(intersect(S, Pb))
    ora_rows[[length(ora_rows) + 1L]] <- data.frame(cohort = co, selector = sel, K = K,
      collection = unname(set_coll[gs]), gs_name = gs, m = m, k = k, nS = nS,
      p = stats::phyper(k - 1L, m, N - m, nS, lower.tail = FALSE), stringsAsFactors = FALSE)
  }
  for (scope in c("cohort", "global")) {
    D  <- if (scope == "cohort") intg$by_type[[intg$cohort_map[[co]]]] else intg$global
    Db <- intersect(D, B); dnum <- length(Db)
    ov <- intersect(S, Db); k <- length(ov)
    drv_rows[[length(drv_rows) + 1L]] <- data.frame(cohort = co, selector = sel, K = K,
      scope = scope, drivers_bg = dnum, overlap = k, nS = nS, nmir = nmir,
      p = stats::phyper(k - 1L, dnum, N - dnum, nS, lower.tail = FALSE),
      genes = paste(ov, collapse = ";"), stringsAsFactors = FALSE)
  }
}
ora <- do.call(rbind, ora_rows); drv <- do.call(rbind, drv_rows)
ora$q <- stats::p.adjust(ora$p, "BH")                          # global BH within ORA family
drv$q <- stats::ave(drv$p, drv$scope, FUN = function(p) stats::p.adjust(p, "BH"))
utils::write.csv(ora, file.path(OUT, "ora_full.csv"), row.names = FALSE)
utils::write.csv(drv, file.path(OUT, "driver_overlap.csv"), row.names = FALSE)

## ---- selector comparison at each K --------------------------------------------
sels <- c("survfcbf_topk", setdiff(unique(sig$selector), "survfcbf_topk"))
summ <- list()
for (K in sort(unique(sig$K))) {
  o <- ora[ora$K == K, ]; dvc <- drv[drv$K == K & drv$scope == "cohort", ]
  row <- lapply(sels, function(sl) {
    oo <- o[o$selector == sl, ]; dd <- dvc[dvc$selector == sl, ]
    enr_per_coh <- tapply(oo$q < Q_SIG, oo$cohort, sum)
    touched     <- tapply(oo$k >= 1L, oo$cohort, sum)          # pathways touched (diversity)
    data.frame(K = K, selector = sl,
      enriched_paths = round(mean(enr_per_coh), 2),
      coh_any_path   = sum(enr_per_coh >= 1),
      paths_touched  = round(mean(touched), 1),
      driver_overlap = round(mean(dd$overlap), 2),
      coh_any_driver = sum(dd$overlap >= 1),
      coh_sig_driver = sum(dd$q < Q_SIG),
      mean_mirna     = round(mean(dd$nmir), 2), stringsAsFactors = FALSE)
  })
  summ[[as.character(K)]] <- do.call(rbind, row)
}
summ <- do.call(rbind, summ); rownames(summ) <- NULL
utils::write.csv(summ, file.path(OUT, "enrichment_summary.csv"), row.names = FALSE)

cat(sprintf("\n=== enrichment summary (K=%d primary; driver overlap = cohort-specific intOGen) ===\n", K_FOCUS))
print(summ[summ$K == K_FOCUS, ], row.names = FALSE)
cat("\n=== K-sensitivity (survfcbf_topk) ===\n")
print(summ[summ$selector == "survfcbf_topk", ], row.names = FALSE)

## ---- figure: driver overlap + enriched pathways by selector (K_FOCUS) ----------
sK <- summ[summ$K == K_FOCUS, ]; sK <- sK[order(-sK$driver_overlap), ]
grDevices::png(file.path(OUT, "enrichment_barplots.png"), width = 1200, height = 520, res = 130)
op <- graphics::par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
cols <- ifelse(sK$selector == "survfcbf_topk", "firebrick", "grey60")
graphics::barplot(sK$driver_overlap, names.arg = sK$selector, las = 2, col = cols,
                  main = sprintf("intOGen driver overlap (K=%d)", K_FOCUS), ylab = "mean drivers in signature")
graphics::barplot(sK$enriched_paths, names.arg = sK$selector, las = 2, col = cols,
                  main = sprintf("Enriched H+C6 sets (BH q<%.2f)", Q_SIG), ylab = "mean enriched pathways")
graphics::par(op); grDevices::dev.off()

saveRDS(list(ora = ora, drv = drv, summ = summ, loss = loss), file.path(OUT, "enrichment.rds"))
cat(sprintf("\ntables + figure -> %s\n", OUT))
cat("ENRICH_DONE\n")
