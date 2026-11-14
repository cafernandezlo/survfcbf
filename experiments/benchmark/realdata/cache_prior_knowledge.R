# Freeze the enrichment prior knowledge as local .rds snapshots (the analysis pipeline
# must not hit a live external server; privacy-safe + reproducible).
# Run ONCE, needs internet for the intOGen download; enrich_signatures.R then reads only
# the cached .rds. Idempotent: skips a cache that already exists.
bench <- Sys.getenv("BENCH_DIR", unset = "")
if (nzchar(bench) && dir.exists(bench)) setwd(bench)
PK <- "realdata/prior_knowledge"
dir.create(PK, recursive = TRUE, showWarnings = FALSE)

## MSigDB Hallmarks (H) + oncogenic (C6) -- msigdbr bundles the sets, no external call.
msig_rds <- file.path(PK, "msigdb_H_C6.rds")
if (!file.exists(msig_rds)) {
  suppressMessages(library(msigdbr))
  keep <- c("gs_name", "gene_symbol", "ensembl_gene", "ncbi_gene")
  grab <- function(coll) {
    df <- as.data.frame(msigdbr(species = "Homo sapiens", collection = coll))[, keep]
    df$collection <- coll; df
  }
  msig <- rbind(grab("H"), grab("C6"))
  saveRDS(list(sets = msig, meta = list(
    msigdbr_version = as.character(utils::packageVersion("msigdbr")),
    date = as.character(Sys.Date()))), msig_rds)
  cat("cached MSigDB H+C6 ->", msig_rds, "\n")
}

## intOGen 2020 Compendium of Cancer Genes -- open download, cached once. Builds a global
## driver set + per-cancer-type sets + the TCGA(Bommert)->intOGen type map for the
## cohort-specific overlap test.
intg_rds <- file.path(PK, "intogen_drivers.rds")
if (!file.exists(intg_rds)) {
  url <- "https://www.intogen.org/download?file=IntOGen-Drivers-20200201.zip"
  zp  <- file.path(PK, "intogen_2020.zip")
  if (!file.exists(zp)) utils::download.file(url, zp, mode = "wb", quiet = TRUE)
  tsv <- grep("Compendium_Cancer_Genes.tsv$", utils::unzip(zp, list = TRUE)$Name, value = TRUE)[1]
  utils::unzip(zp, files = tsv, exdir = PK)
  d <- utils::read.delim(file.path(PK, tsv), stringsAsFactors = FALSE)
  d$SYMBOL <- trimws(d$SYMBOL)
  cohort_map <- c(BLCA = "BLCA", BRCA = "BRCA", HNSC = "HNSC", KIRC = "RCCC", LGG = "LGG",
                  LUAD = "LUAD", LUSC = "LUSC", OV = "OV", PAAD = "PAAD", SKCM = "CM", STAD = "ST")
  saveRDS(list(global = sort(unique(d$SYMBOL)),
               by_type = lapply(split(d$SYMBOL, d$CANCER_TYPE), function(x) sort(unique(x))),
               cohort_map = cohort_map,
               source = "intOGen 2020-02-02 (IntOGen-Drivers-20200213)",
               date = as.character(Sys.Date())), intg_rds)
  cat("cached intOGen drivers ->", intg_rds, "\n")
}

## HGNC protein-coding gene set -- open download, cached once. Drives the coding-only
## enrichment sub-analysis (extract_signatures.R CODING_ONLY mode restricts to these
## Ensembl ids; ~71% of the Bommert rna features are protein-coding).
pc_rds <- file.path(PK, "protein_coding.rds")
if (!file.exists(pc_rds)) {
  url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_groups/protein-coding_gene.txt"
  tsv <- file.path(PK, "hgnc_protein_coding.txt")
  if (!file.exists(tsv)) utils::download.file(url, tsv, mode = "wb", quiet = TRUE)
  h <- utils::read.delim(tsv, stringsAsFactors = FALSE, quote = "")
  saveRDS(list(ensembl = unique(h$ensembl_gene_id[nzchar(h$ensembl_gene_id)]),
               symbol = unique(h$symbol), source = "HGNC protein-coding_gene.txt",
               date = as.character(Sys.Date())), pc_rds)
  cat("cached HGNC protein-coding ->", pc_rds, "\n")
}
cat("prior knowledge ready in", PK, "\n")
