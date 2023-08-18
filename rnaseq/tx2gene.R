#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --mem=16G
#SBATCH --job-name=tx2gene
#SBATCH --output=slurm-tx2gene-%j.out

# DESCRIPTION ------------------------------------------------------------------
# This script will create a transcript-("tx")-to-gene TSV file from a GTF
# Version 2023-08-16 by Jelmer Poelstra

# With the following Conda environment, the script wouldn't need to install anything:
# /fs/ess/PAS0471/jelmer/conda/r-rnaseq

#TODO - Be more flexible about identifier selection, this will probably fail for some GTF files
#TODO - Accept GFF as well

# SET-UP -----------------------------------------------------------------------
# Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(ape)) install.packages("ape", repos = rep, lib = lib)
  if (!require(tidyverse)) install.packages("tidyverse", repos = rep, lib = lib)
} )
library(argparse)
library(ape)
library(tidyverse)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--gtf",
                    type = "character", required = TRUE, default = NULL,
                    help = "Input GTF file (REQUIRED)")
parser$add_argument("-o", "--tx2gene",
                    type = "character", required = TRUE, default = NULL,
                    help = "Output transcript-to-gene file")
args <- parser$parse_args()

gtf_file <- args$gtf            # gtf_file <- "data/ref/GCF_000001405.40_ed.gtf"
tx2gene_file <- args$tx2gene

# Define output files
outdir <- dirname(tx2gene_file)
gene_tab_file <- file.path(outdir, "gene_table.tsv")

# Report
message("\n# Starting script tx2gene.R")
Sys.time()
message()
message("Input GTF file:                ", gtf_file)
message("Output tx2gene file:           ", tx2gene_file)
message("Output gene-table file:        ", gene_tab_file)
message("======================================================================")
message()

# CREATE THE TX2GENE TABLE -----------------------------------------------------
# Create the output dir
dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Read and process the GTF files
gtf_df <- read.gff(gtf_file, GFF3 = FALSE) |> 
  filter(feature == "transcript") |> 
  mutate(
    gene_id = sub(".*gene_id \"([^\"]+)\".*", "\\1", attributes),
    transcript_id = sub(".*transcript_id \"([^\"]+)\".*", "\\1", attributes),
    db_xref = sub(".*db_xref \"([^\"]+)\".*", "\\1", attributes),
    description = sub(".*product \"([^\"]+)\".*", "\\1", attributes),
  ) |>
  select(gene_id, transcript_id, db_xref, description)

# Create a transcript-to-gene df (1 row per transcript)
tx2gene <- gtf_df |> select(transcript_id, gene_id)

# Create a df with per-gene info (1 row per gene, no transcripts)
gene_tab <- gtf_df |>
  select(gene_id, db_xref, description) |>
  distinct(gene_id, .keep_all = TRUE)

# Write output files
write_tsv(tx2gene, tx2gene_file, col_names = FALSE)
write_tsv(gene_tab, gene_tab_file, col_names = FALSE)

# Report
message("Nr of rows in the tx2gene table (= nr transcripts):    ", nrow(tx2gene))
message("Nr of rows in the gene-info table (= nr genes):        ", nrow(gene_tab))

# WRAP UP ----------------------------------------------------------------------
# List output
message("\n# Listing the output file(s):")
system(paste("ls -lh", tx2gene_file))
system(paste("ls -lh", gene_tab_file))
message("\n# Done with script tx2gene.R")
Sys.time()
message()
sessionInfo()

# Alternative way - From https://www.biostars.org/p/9485199/
#? HOWEVER, this removed three transcripts with bad stop codons/CDS, while we want to keep all
#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
#k <- keys(txdb, keytype = "TXNAME")
#tx2gene <- select(txdb, k, "GENEID", "TXNAME")
