#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=90
#SBATCH --mem=40G
#SBATCH --cpus-per-task=10
#SBATCH --job-name=strandcheckr
#SBATCH --output=slurm-strandcheckr-%j.out


# DESCRIPTION ------------------------------------------------------------------
# This script will try to remove DNA contamination from RNAseq BAM files,
# using the 'strandCheckR' R package. See:
# - https://uofabioinformaticshub.github.io/strandCheckR/
# - https://bioconductor.org/packages/release/bioc/vignettes/strandCheckR/inst/doc/strandCheckR.html


# SETUP ------------------------------------------------------------------------
# Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse, quietly = TRUE)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman, quietly = TRUE)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(BiocManager, quietly = TRUE)) install.packages("BiocManager", repos = rep, lib = lib)
  if (!require(strandCheckR, quietly = TRUE)) install.packages("tidyverse", repos = rep, lib = lib)
  if (!require(ggplot2, quietly = TRUE)) BiocManager::install("DESeq2")
  packages <- c("argparse", "DESeq2", "tidyverse")
  pacman::p_load(char = packages, install = TRUE)
} )

# Parse options
parser <- ArgumentParser()
parser$add_argument("-i", "--bam",
                    type = "character",
                    required = FALSE,
                    default = "results/nfc_rnaseq",
                    help = "Input dir with nf-core RNAseq workflow results")
parser$add_argument("-o", "--outdir",
                    type = "character",
                    required = FALSE,
                    default = "dds.rds",
                    help = "Output filename (can be in another dir)")
args <- parser$parse_args()
bam <- args$bam
outdir <- args$outdir

# Define the output files
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
file_id <- sub(".markdup.sorted.bam", "", basename(bam))
bam_out <- file.path(outdir, paste0(file_id, ".bam"))
stat_out <- file.path(outdir, paste0(file_id, "_stats.txt"))
fig_out <- file.path(outdir, paste0(file_id, "_hist.png"))

# Report
message("# Starting script strandcheckr.R")
Sys.time()
message()
message("# Input BAM file:                   ", bam)
message("# Output dir:                       ", outdir)
message("======================================================================")
message()


# RUN --------------------------------------------------------------------------
# Filter the BAM file
filter_df <- filterDNA(
  file = bam,
  destination = bam_out,
  statFile = stat_out,
  paired = TRUE,
  winWidth = 500,
  threshold = 0.7,
  getWin = TRUE
)

# Create a plot
filter_df$File <- basename(as.character(filter_df$File))
p <- plotHist(filter_df,
              groupBy = "File", normalizeBy = "File", scales = "free_y")
ggsave(fig_out, p, width = 8, height = 6.5, dpi = "retina")


# WRAP UP ----------------------------------------------------------------------
# Show summary stats
message("\n# Printing filtering summary stats:")
system(paste("cat", stat_out))

# List output files
message("\n# Listing the output file(s):")
system(paste("ls -lh", bam_out))
system(paste("ls -lh", stat_out))
system(paste("ls -lh", fig_out))

message("\n# Done with script strandcheckr.R")
Sys.time()
message()
sessionInfo()
