#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --mem=16G
#SBATCH --job-name=tx2gene
#SBATCH --output=slurm-tx2gene-%j.out

# DESCRIPTION ------------------------------------------------------------------
# This script will create a transcript-("tx")-to-gene TSV file from a GTF
# Version 2023-12-165 by Jelmer Poelstra

# With the following Conda environment, the script wouldn't need to install anything:
# /fs/ess/PAS0471/jelmer/conda/r-rnaseq

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

# Parse command-line arguments passed to the script
parser <- ArgumentParser()
parser$add_argument("-i", "--gtf",
                    type = "character", required = TRUE, default = NULL,
                    help = "Input GTF file (REQUIRED)")
parser$add_argument("-o", "--tx2gene",
                    type = "character", required = TRUE, default = NULL,
                    help = "Output transcript-to-gene file")
args <- parser$parse_args()
gtf_file <- args$gtf
tx2gene_file <- args$tx2gene

# Define output files
outdir <- dirname(tx2gene_file)
gene_tab_file <- file.path(outdir, "gene_table.tsv")

# Report
message("\n# Starting script tx2gene.R")
Sys.time()
message()
message("# Input GTF file:                ", gtf_file)
message("# Output tx2gene file:           ", tx2gene_file)
message("# Output gene-table file:        ", gene_tab_file)
message("======================================================================")
message()

# CREATE THE TX2GENE TABLE -----------------------------------------------------
# Create the output dir
dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Read the GTF file and select transcripts
gtf_df_raw <- read.gff(gtf_file, GFF3 = FALSE) |> 
  filter(feature == "transcript")

# Extract transcript info
gtf_df <- gtf_df_raw |>
  mutate(
    gene_id = str_match(attributes, ".*gene_id \"([^\"]+)\".*")[, 2],
    transcript_id = str_match(attributes, ".*transcript_id \"([^\"]+)\".*")[, 2],
    db_xref = str_match(attributes, ".*db_xref \"([^\"]+)\".*")[, 2],
    gene_name = str_match(attributes, ".*gene_name \"([^\"]+)\".*")[, 2],
    biotype = str_match(attributes, ".*gene_biotype \"([^\"]+)\".*")[, 2],
    description = str_match(attributes, ".*product \"([^\"]+)\".*")[, 2]
  ) |>
  select(gene_id, transcript_id, biotype, gene_name, db_xref, description)

# Stop if there are no gene or transcript IDs
stopifnot("No gene IDs found" = !all(is.na(gtf_df$gene_id)))
stopifnot("No transcript IDs found" = !all(is.na(gtf_df$transcript_id)))

# Remove columns that only have NA
if (all(is.na(gtf_df$biotype))) gtf_df$biotype <- NULL
if (all(is.na(gtf_df$db_xref))) gtf_df$db_xref <- NULL
if (all(is.na(gtf_df$description))) gtf_df$description <- NULL

# Create a transcript-to-gene df (1 row per transcript)
tx2gene <- gtf_df |> select(transcript_id, gene_id)
message("First few rows of the tx2gene dataframe:")
print(head(tx2gene))

# Create a df with per-gene info (1 row per gene, no transcripts)
gene_tab <- gtf_df |>
  select(any_of(c("gene_id", "biotype", "gene_name", "db_xref", "description"))) |>
  distinct(gene_id, .keep_all = TRUE)
message("\nFirst few rows of the gene info dataframe:")
print(head(gene_tab))

# Write output files
write_tsv(tx2gene, tx2gene_file, col_names = FALSE)
write_tsv(gene_tab, gene_tab_file, col_names = FALSE)

# Report
message("\n# Nr of rows in the tx2gene table (= nr transcripts):    ", nrow(tx2gene))
message("# Nr of rows in the gene-info table (= nr genes):        ", nrow(gene_tab))

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
