#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --job-name=ggtree
#SBATCH --output=slurm-ggtree-%j.out

#? From an input tree file, this script will plot the tree with ggtree
#? (All tree file formats should be supported)

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
parser$add_argument("-i", "--tree",
                    type = "character", required = TRUE,
                    help = "Input tree file (REQUIRED)")
parser$add_argument("-o", "--figure",
                    type = "character", required = FALSE, default = NULL,
                    help = "Output figure file in '.png' format
                            [default: same dir and name as tree file, but .png extension]")
args <- parser$parse_args()
tree_file <- args$tree
figure_file <- args$figure

# Define the output file name, if needed
if (is.null(figure_file)) {
  figure_file <- file.path(
    dirname(tree_file),
    paste0(tools::file_path_sans_ext(basename(tree_file)), ".png")
  )
}

# Report
message("\n# Starting script ggtree.R")
message("# Input tree file:          ", tree_file)
message("# Output figure file:       ", figure_file, "\n")


# PLOT THE TREE ----------------------------------------------------------------
# Read the tree file and prep the tree
tree <- read.tree(tree_file)
nseqs <- length(tree$tip.label)
total_edge_len <- sum(tree$edge.length)

#tree <- ape::root(tree, outgroup = "AG21-0056") #TODO - Add 'root' option

# Make the plot
p <- ggtree(tree, layout = "rectangular") +
        geom_tiplab(size = 2.5, color = "grey50") +
        geom_treescale(x = 0, y = nseqs - 1, color = "grey50") +
        geom_rootedge(rootedge = sum(tree$edge.length) / 50) +
        theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm")) +
        coord_cartesian(clip = "off")

# Save the plot to file
ggsave(figure_file, p, width = 6, height = 6, dpi = "retina")


# WRAP UP ----------------------------------------------------------------------
message("# Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n# Done with script ggtree.R")
Sys.time()
