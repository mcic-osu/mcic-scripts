#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --mem=170G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=tximport
#SBATCH --output=slurm-tximport-%j.out


# DESCRIPTION ------------------------------------------------------------------
# This script will import Kallisto files into R and create a tximport object,
# and if metadata is provided, also a DEseq2 object
#TODO - Also accept Salmon counts


# SET-UP -----------------------------------------------------------------------
# Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(BiocManager)) install.packages("BiocManager", repos = rep, lib = lib)
  if (!require(tximport)) BiocManager::install("tximport")
  if (!require(rhdf5)) BiocManager::install("rhdf5")
  if (!require(DESeq2)) BiocManager::install("DESeq2")
} )
packages <- c("argparse", "tximport", "rhdf5", "DESeq2", "tidyverse")
pacman::p_load(char = packages, install = TRUE)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--indir",
                    type = "character", required = TRUE, default = NULL,
                    help = "Input dir with Kallisto output files (REQUIRED)")
parser$add_argument("-o", "--outdir",
                    type = "character", required = TRUE, default = NULL,
                    help = "Output directory")
parser$add_argument("--tx2gene",
                    type = "character", required = TRUE, default = NULL,
                    help = "Transcript-to-gene map (TSV)")
parser$add_argument("--meta",
                    type = "character", required = FALSE, default = NULL,
                    help = "Metadata file (TSV), needed to create a DESeq object")
parser$add_argument("--id_col",
                    type = "integer",
                    required = FALSE,
                    default = 1,
                    help = "In metadata file, column number with sample IDs")
args <- parser$parse_args()

indir <- args$indir
outdir <- args$outdir
tx2gene_file <- args$tx2gene
meta_file <- args$meta
sample_id_column <- args$id_col

# Define output files
txi_out <- file.path(outdir, "tximport_object.rds")
dds_out <- file.path(outdir, "deseq_object.rds")

# Report
message("\n# Starting script tximport.R")
Sys.time()
message()
message("Input dir:                     ", indir)
message("Output dir:                    ", outdir)
message("Transcript-to-gene map file:   ", tx2gene_file)
message("======================================================================")
message()


# PROCESS KALLISTO COUNTS ------------------------------------------------------
# Create output dir
dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Read transcript-to-gene map
tx2gene <- read_tsv(tx2gene_file,
                     col_names = c("GENEID", "TXNAME"),
                     col_types = "cc") |>
  dplyr::select(TXNAME, GENEID)

# Import Kallisto transcript counts --
# create gene-level count estimates normalized by library size and transcript length
# See https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
kallisto_files <- list.files(indir, pattern = "abundance.h5$", recursive = TRUE,
                             full.names = TRUE)

txi <- tximport(kallisto_files,
                type = "kallisto",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

# Save tximport object
saveRDS(txi, txi_out)


# CREATE DESEQ OBJECT ----------------------------------------------------------
if (!is.null(meta_file)) {
    # Read metadata
    meta <- read.delim(meta_file) |>
        arrange(.data[[sample_id_column]])
    rownames(meta) <- meta[[sample_id_column]]

    # Name Kallisto files according to metadata
    names(kallisto_files) <- meta[[sample_id_column]]

    # Check that sample names are the same, and that samples are in the same order
    stopifnot(all(rownames(meta) == colnames(txi$counts)))
    message("\n# Sample names:")
    print(rownames(meta))
    message("\n# Dimensions of count matrix:")
    dim(txi$counts)

    # Create DESeq object
    dds <- DESeqDataSetFromTximport(txi, meta, ~1)

    # Save DESeq object
    saveRDS(dds, dds_out)
}


# WRAP UP ----------------------------------------------------------------------
# List output
message("\n# Listing the output file(s):")
system(paste("ls -lh", txi_out))
if (!is.null(meta_file)) system(paste("ls -lh", dds_out))

message("\n# Done with script tximport.R")
Sys.time()
message()
sessionInfo()
