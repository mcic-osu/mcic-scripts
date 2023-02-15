#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=10
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=make_dds
#SBATCH --output=slurm-make_dds-%j.out


# DESCRIPTION ------------------------------------------------------------------
# This script will create a DESeq object from one of the RDS files that the
# nf-core RNAseq pipeline produces.


# SETUP ------------------------------------------------------------------------
# Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse, quietly = TRUE)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman, quietly = TRUE)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager", repos = rep, lib = lib)
  if (!require(tidyverse, quietly = TRUE)) install.packages("tidyverse", repos = rep, lib = lib)
  if (!require(DESeq2, quietly = TRUE)) BiocManager::install("DESeq2")
  packages <- c("argparse", "rhdf5", "DESeq2", "tidyverse")
  pacman::p_load(char = packages, install = TRUE)
} )

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--indir",
                    type = "character",
                    required = FALSE,
                    default = "results/nfc_rnaseq",
                    help = "Input dir with nf-core RNAseq workflow results")
parser$add_argument("-o", "--outfile",
                    type = "character",
                    required = FALSE,
                    default = "dds.rds",
                    help = "Output filename (can be in another dir)")
parser$add_argument("--meta",
                    type = "character",
                    required = FALSE,
                    default = NULL,
                    help = "Metadata file (TSV)")
parser$add_argument("--id_col",
                    type = "integer",
                    required = FALSE,
                    default = 1,
                    help = "In metadata file, column number with sample IDs")
args <- parser$parse_args()

indir <- args$indir
outfile <- args$outdir
meta_file <- args$meta
sample_id_column <- args$id_col

# Define input files
infile_part <- "star_salmon/salmon.merged.gene_counts_length_scaled.rds"
infile <- file.path(indir, infile_rds_part)

# Output files
outdir <- dirname(outfile)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Report
message("Starting script nfc-rnaseq_make-deseq.R")
Sys.time()
message()
message("Input dir:                        ", indir)
message("Input RDS file:                   ", infile_rds)
message("Output RDS file (DESeq object):   ", outfile)
if (!is.null(meta_file)) message("Input metadata file:              ", meta_file)
message("======================================================================")
message()


# CREATE THE DESEQ OBJECT ------------------------------------------------------
# Load counts
count_obj <- readRDS(infile_rds)
count_mat <- round(assays(count_obj)$counts)

# Metadata
if (!is.null(meta_file)) {
  meta_df <- read.delim(metadata_file, header = TRUE, sep = "\t")
  
  # Sample IDs as rownames
  rownames(meta_df) <- meta_df[[id_col]]
  
  # Sort by sample ID
  meta_df <- meta_df[order(meta_df[[1]]), ]
  
  # Make columns factors (except the 1st one, which should be sample IDs)
  cols <- colnames(meta_df[2:ncol(meta_df)])
  meta_df[cols] <- lapply(meta_df[cols], as.factor)
  
  # Check that sample names are the same, and that samples are in the same order
  stopifnot(all(rownames(meta) == colnames(count_mat)))
  
  # Report
  message("\n# Sample names:")
  print(rownames(meta))
  message("\n# Dimensions of count matrix:")
  dim(count_mat)
  
} else {
  meta_df <- counts@colData
}

# Create DESeq object (For now, set design to dummy `~1` => intercept only)
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = meta_df,
                              design = ~1)

# Save the output file
saveRDS(dds, outfile)


# WRAP UP ----------------------------------------------------------------------
# List output
message("\n# Listing the output file(s):")
system(paste("ls -lh", outfile))

message("\n# Done with script nfc-rnaseq_make-deseq.R")
Sys.time()
message()
sessionInfo()

# Docs:
#  - https://nf-co.re/rnaseq/3.10.1/output#salmon
#  - https://github.com/nf-core/rnaseq/issues/499
#  - https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html

# Alternatively, could load RData from the workflow's DESeq analysis (these are the same)
# infile_deseq <- file.path(NFCORE_RNASEQ_OUTDIR, "star_salmon/deseq2_qc/deseq2.dds.RData")
# load(infile_deseq)
