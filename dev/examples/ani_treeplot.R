#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --job-name=ani_treeplot
#SBATCH --output=slurm-ani_treeplot-%j.out

#? From an input tree file, this script will plot a time-tree with ggtree
#? (All tree file formats should be supported)

#? Load the Conda environment as follows to run this script without it needing to install R packages:
#? module load miniconda3 && source activate /fs/ess/PAS0471/jelmer/conda/r_tree

# SET-UP -----------------------------------------------------------------------
# Load packages
packages <- c("argparse", "ggtree", "ggplot2", "ape")
pacman::p_load(char = packages, install = TRUE)

# Variables: parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--ani_csv",
                    type = "character", required = TRUE,
                    help = "Input CSV file with ANI or equivalent similarity scores (REQUIRED)")
parser$add_argument("-o", "--figure",
                    type = "character", required = FALSE, default = NULL,
                    help = "Output figure file in '.png' format
                            [default: same dir and name as ANI file, but .png extension]")
args <- parser$parse_args()

# Save args from args list into separate variables
ani_file <- args$ani_csv
figure_file <- args$figure

# Define the output file name, if needed
if (is.null(figure_file)) {
  figure_file <- file.path(
    dirname(ani_file),
    paste0(tools::file_path_sans_ext(basename(ani_file)), ".png")
  )
}

# Report
message("\n# Starting script ani_treeplot.png")
message("# Input file:             ", ani_file)
message("# Output figure file:     ", figure_file, "\n")


# READ AND PROCESS THE INPUT FILES ---------------------------------------------
# Read the input files
ani <- read.csv(ani_file)
colnames(ani) <- sub(".fna.sig", "", colnames(ani))

# Get a distance matrix
dist_mat <- as.matrix(1 - ani)
rownames(dist_mat) <- colnames(dist_mat)
dists <- as.dist(dist_mat)

# Build a tree (https://rpubs.com/WalshJake75/674724)
tree <- phangorn::upgma(dists)

# PLOT THE TREE ----------------------------------------------------------------
# Build a tree with ggtree
p <- ggtree(tree, layout = "circular") +
  geom_tiplab(size = 2, color = "grey50") +
  geom_rootedge(rootedge = sum(tree$edge.length) / 100) +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))

ggsave(figure_file, p, width = 8, height = 8, dpi = "retina")


# WRAP UP ----------------------------------------------------------------------
message("\n# Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n# Done with script ani_treeplot.R")
Sys.time()

# NOTES ------------------------------------------------------------------------
# Find node numbers (https://bioconnector.github.io/workshops/r-ggtree.html):
# ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

# Could also build tree with (but this resulted in negative branch lengths):
#tree <- ape::nj(dists)

# And a useful quick ape plot of the tree:
#plot(tree, "unrooted")
