#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=tree_build
#SBATCH --output=slurm-tree_build-%j.out

#? This script will build a phylogenetic tree from a set of ASVs

#? Load the Conda environment as follows to run this script directly using sbatch:
#? module load miniconda3/4.12.0-py39 && source activate /fs/ess/PAS0471/jelmer/conda/r-metabar

# SET-UP -----------------------------------------------------------------------
# Packages
packages <- c("BiocManager", "dada2", "DECIPHER", "phangorn")

# Get nr of cores
n_threads <- as.integer(system("echo $SLURM_CPUS_PER_TASK", intern = TRUE))

# Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)
parser <- ArgumentParser()
parser$add_argument("-i", "--infile",
                    type = "character", required = TRUE,
                    help = "Input file (sequence table RDS or FASTA file) (REQUIRED)")
parser$add_argument("-o", "--tree",
                    type = "character", required = TRUE,
                    help = "Output file (tree RDS file) (REQUIRED)")
args <- parser$parse_args()

infile <- args$infile
tree_rds <- args$tree

# Infer input format (FASTA or not)
if (grepl(".fas?t?a?$", infile)) fasta_format <- TRUE else fasta_format <- FALSE

# Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(char = packages)

# Create output dir if needed
outdir <- dirname(tree_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Report
message()
message("# ====================================================================")
message("#               STARTING SCRIPT TREE_BUILD.R")
message("# ====================================================================")
Sys.time()
message()
message("# Input file:                         ", infile)
message("# Input file is in FASTA format:      ", fasta_format)
message("# Output file (tree RDS file):        ", tree_rds)
message("# Number of cores:                    ", n_threads)
message("# ====================================================================")
message()


# LOAD DATA --------------------------------------------------------------------
if (fasta_format == FALSE) {
  # Assume input file is a sequence table from dada
  seqtab <- readRDS(infile)
  seqs <- getSequences(seqtab)
  names(seqs) <- seqs
} else {
  # Assume input file is a FASTA file
  seqs <- readDNAStringSet(infile)
}


# BUILD THE TREE ---------------------------------------------------------------
# See https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html

# 1 - Align
message("\n# Aligning sequences...")
alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor = NA,
                       iterations = 5,
                       refinements = 5,
                       processors = n_threads)

# 2 - Compute distances
message("\n# Computing pairwise distances from ASVs...")
alignment_mat <- phyDat(as(alignment, "matrix"), type = "DNA")
dist_mat <- dist.ml(alignment_mat)

# 3 - Build neighbor-joining tree and compute its likelihood
message("\n# Building a tree...")
nj_tree <- NJ(dist_mat) # Build NJ tree
fit <- pml(nj_tree, data = alignment_mat) # Compute likelihood
fit_gtr <- update(fit, k = 4, inv = 0.2) # Update to GTR model

# 4 - Compute likelihood
message("\n# Optimizing the tree...")
fit_gtr <- optim.pml(fit_gtr,
                     model = "GTR",
                     optInv = TRUE,
                     optGamma = TRUE,
                     rearrangement = "stochastic",
                     control = pml.control(trace = 0))

# 5 - Root the tree
message("\n# Rooting the tree...")
rooted_tree <- root_by_longest_edge(fit_gtr)


# WRAP UP ----------------------------------------------------------------------
# Save tree RDS file
saveRDS(rooted_tree, tree_rds)

# Report
message("\n# Listing output file:")
system(paste("ls -lh", tree_rds))
message()
message("# Done with script tree.R")
Sys.time()
message()
sessionInfo()
message()
