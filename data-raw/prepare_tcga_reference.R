# ============================================================================
# Prepare the TCGA BRCA empirical-marginal reference for the simulator.
#
# Reads the TCGAbiolinks STAR-Counts SummarizedExperiment, extracts the
# tpm_unstrand assay, applies log2(x + 1), restricts to expressed protein-coding
# genes (collapsing duplicate gene symbols), and writes a plain genes-by-samples
# matrix (no Bioconductor classes) consumable by build_empirical_reference().
#
# Bioconductor (SummarizedExperiment) is needed ONLY here, not by the package.
# The TCGA inputs and the produced reference are sensitive and local-only: they
# live outside the repo and never leave the machine. Run once.
# ============================================================================

suppressMessages(library(SummarizedExperiment))

tcga_dir <- Sys.getenv("TCGA_DIR", "")
if (!nzchar(tcga_dir))
  stop("set TCGA_DIR to the local folder holding tcga_brca_se.rds", call. = FALSE)
if (!dir.exists(tcga_dir)) stop("TCGA_DIR does not exist: ", tcga_dir, call. = FALSE)

in_file  <- file.path(tcga_dir, "tcga_brca_se.rds")
out_file <- file.path(tcga_dir, "tcga_brca_log2tpm_reference.rds")

message("Reading: ", in_file)
se <- readRDS(in_file)

tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")
rd  <- as.data.frame(SummarizedExperiment::rowData(se))
rm(se); gc()

mat <- log2(tpm + 1)
rm(tpm)

## protein-coding genes only, labelled by HGNC symbol
keep <- rd$gene_type == "protein_coding"
mat <- mat[keep, , drop = FALSE]
rownames(mat) <- rd$gene_name[keep]

## drop unexpressed or invariant genes
expressed <- rowMeans(mat > 0) >= 0.5
row_var   <- apply(mat, 1L, stats::var)
mat <- mat[expressed & is.finite(row_var) & row_var > 0, , drop = FALSE]

## collapse duplicate gene symbols, keeping the most-expressed copy
if (anyDuplicated(rownames(mat))) {
  ord <- order(rowMeans(mat), decreasing = TRUE)
  mat <- mat[ord, , drop = FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
}

message(sprintf("BRCA log2(TPM+1) reference: %d genes x %d samples", nrow(mat), ncol(mat)))
saveRDS(mat, out_file)
message("Saved: ", out_file)
