#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-tree-%j.out


# SET-UP -----------------------------------------------------------------------
## Report
cat("## Starting script tree.R\n")
Sys.time()
message()

## Parse command-line arguments
if(!"argparse" %in% installed.packages()) install.packages("argparse")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--seqtab",
                    type = "character", required = TRUE,
                    help = "Input file (sequence table RDS) (REQUIRED)")
parser$add_argument("-o", "--tree",
                    type = "character", required = TRUE,
                    help = "Output file (tree RDS file) (REQUIRED)")
args <- parser$parse_args()

seqtab_rds <- args$seqtab
tree_rds <- args$tree

## Other variables
n_cores <- as.integer(system("echo $SLURM_CPUS_PER_TASK", intern = TRUE))

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("BiocManager", "dada2", "DECIPHER", "phangorn")
pacman::p_load(char = packages)

## Report
message("## Input file (sequence table RDS):    ", seqtab_rds)
message("## Output file (tree RDS file):        ", tree_rds)
message()
message("## Number of cores:                    ", n_cores)
message()

## Create output dir if needed
outdir <- dirname(tree_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Load the input data
seqtab <- readRDS(seqtab_rds)
seqs <- getSequences(seqtab)
names(seqs) <- seqs


# BUILD THE TREE ---------------------------------------------------------------
# See https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html

## 1 - Align
message("\n## Aligning sequences...")
alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor = NA,
                       iterations = 5,
                       refinements = 5,
                       processors = n_cores)

## 2 - Compute distances
message("\n## Computing pairwise distances from ASVs...")
alignment_mat <- phyDat(as(alignment, "matrix"), type = "DNA")
dist_mat <- dist.ml(alignment_mat)

## 3 - Build neighbor-joining tree and compute its likelihood
message("\n## Building a tree...")
nj_tree <- NJ(dist_mat) # Build NJ tree
fit <- pml(nj_tree, data = alignment_mat) # Compute likelihood
fit_gtr <- update(fit, k = 4, inv = 0.2) # Update to GTR model

## 4 - Compute likelihood
message("\n## Optimizing the tree...")
fit_gtr <- optim.pml(fit_gtr,
                     model = "GTR",
                     optInv = TRUE,
                     optGamma = TRUE,
                     rearrangement = "stochastic",
                     control = pml.control(trace = 0))


# WRAP UP ----------------------------------------------------------------------
## Save tree RDS file
saveRDS(fit_gtr, tree_rds)

## Report
message("\n## Listing output file:")
system(paste("ls -lh", tree_rds))
message()
message("## Done with script tree.R")
Sys.time()
message()
