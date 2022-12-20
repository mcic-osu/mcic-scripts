## Packages
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(DESeq2)) {
  BiocManager::install("DESeq2")
  library(DESeq2)
}

## Constants
NFCORE_RNASEQ_OUTDIR <- "results/nfc_rnaseq/"

## Define input files
infile_rds <- file.path(NFCORE_RNASEQ_OUTDIR,
                        "star_salmon/salmon.merged.gene_counts_length_scaled.rds")

## Alternatively, could load RData from the workflow's DESeq analysis (these are the same)
#infile_deseq <- file.path(NFCORE_RNASEQ_OUTDIR,
#                          "star_salmon/deseq2_qc/deseq2.dds.RData")
#load(infile_deseq)

## Load counts
counts <- readRDS(infile_rds)

## Create DESeq object.
## For now, set design to dummy `~1` (intercept only)
dds <- DESeqDataSetFromMatrix(countData = round(assays(counts)$counts),
                              colData = counts@colData,
                              design = ~1)  

## For next steps,
## see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

