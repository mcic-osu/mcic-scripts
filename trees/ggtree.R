#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --job-name=ggtree
#SBATCH --output=slurm-ggtree-%j.out

#? Load the Conda environment as follows to run this script without it needing to install R packages:
#? module load miniconda3 && source activate /fs/ess/PAS0471/jelmer/conda/r_tree
#(micromamba create -y -p /fs/ess/PAS0471/jelmer/conda/r_tree -c bioconda r-argparse r-pacman bioconductor-ggtree r-biocmanager r-ape)

# SET-UP -----------------------------------------------------------------------
packages <- c("ape", "ggtree", "argparse", "ggplot2")
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(BiocManager)) install.packages("BiocManager", repos = rep, lib = lib)
  if (!require(ggtree)) BiocManager::install("ggtree")
  if (!require(ggplot2)) install.packages("ggplot2", repos = rep, lib = lib)
  if (!require(ape)) install.packages("ape", repos = rep, lib = lib)
} )
pacman::p_load(char = packages, install = TRUE)

# Parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--tree_file",
                    type = "character", required = TRUE,
                    help = "Input tree file in '.tree' format (REQUIRED)")
parser$add_argument("-o", "--figure_file",
                    type = "character", required = TRUE,
                    help = "Output figure file in '.png' format (REQUIRED)")
args <- parser$parse_args()
tree_file <- args$tree_file
figure_file <- args$figure_file

# Report
message("\n# Starting script ggtree.R")
message("# Input tree file:          ", tree_file)
message("# Output figure file:       ", figure_file, "\n")


# PLOT THE TREE ----------------------------------------------------------------
# Read the tree file and prep the tree
tree <- read.tree(tree_file)
nseqs <- length(tree$tip.label)

# Make the plot
p <- ggtree(tree, layout = "rectangular") +
        geom_tiplab(size = 2.5) +
        geom_treescale(x = 0, y = nseqs - 1, color = "grey50") +
        geom_rootedge(rootedge = 0.005) +
        theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm")) +
        coord_cartesian(clip = "off")

# Save the plot to a figure
ggsave(figure_file, p)


# WRAP UP ----------------------------------------------------------------------
sessionInfo()
message("\n# Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n# Done with script ggtree.R")
Sys.time()
