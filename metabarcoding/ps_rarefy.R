#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm-ps_rarefy-%j.out

#? NOTE: In the metadata file, sample IDs should be in the 1st column

# SET-UP -----------------------------------------------------------------------
## Report
message()
message("## Starting script ps_rarefy.R")
Sys.time()
message()

## Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--ps_in",
                    type = "character", required = TRUE,
                    help = "Input RDS file with phyloseq object (REQUIRED)")
parser$add_argument("-o", "--ps_out",
                    type = "character", required = TRUE,
                    help = "Output RDS file with phyloseq object (REQUIRED)")
args <- parser$parse_args()

ps_in_file <- args$ps_in
ps_out_file <- args$ps_out

## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "phyloseq")
pacman::p_load(char = packages)

## Create output dir if needed
outdir <- dirname(ps_out_file)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Report
message("## Phyloseq input RDS file:          ", ps_in_file)
message("## Phyloseq output RDS file:         ", ps_out_file)
message()


# LOAD INPUT DATA AND RAREFY ---------------------------------------------------
ps_in <- readRDS(ps_in_file)

ps_out <- rarefy_even_depth(ps_in,
                            sample.size = min(sample_sums(ps_in)),
                            rngseed = 1, replace = FALSE)

message("## Rarefied to a per-sample nr of ASVs of: ", min(sample_sums(ps_in)))
message("## Total nr of ASVs in input object: ", sum(sample_sums(ps_in)))
message("## Total nr of ASVs in output object: ", sum(sample_sums(ps_out)))

# SAVE AND WRAP UP -------------------------------------------------------------
## Save phyloseq object in an RDS file
saveRDS(ps_out, ps_out_file)

## Report
message("\n## Listing output file:")
system(paste("ls -lh", ps_out_file))
message("\n## Done with script ps_rarefy.R")
Sys.time()
message()
