#!/usr/bin/env Rscript

# SET-UP -----------------------------------------------------------------------
# remotes::install_version("dplyr", version = "1.0.5")
# remotes::install_github("YuLab-SMU/tidytree")
# remotes::install_github("YuLab-SMU/ggtree")
if (!"pacman" %in% installed.packages())
  install.packages("pacman", repos = "https://cloud.r-project.org")
if (!"BiocManager" %in% installed.packages())
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!"ggtree" %in% installed.packages()) BiocManager::install("ggtree")
packages <- c("tidyverse", "here", "ape", "ggtree", "argparse")
pacman::p_load(char = packages, install = TRUE, repos = "https://cloud.r-project.org")

## Parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--tree_file",
                    type = "character", required = TRUE,
                    help = "Input tree file in .tree format (REQUIRED)")
parser$add_argument("-o", "--figure_file",
                    type = "character", required = TRUE,
                    help = "Output figure file in .png format (REQUIRED)")
args <- parser$parse_args()

tree_file <- args$tree_file
figure_file <- args$figure_file

## Report
message("\n## Starting script ggtree.R")
message("## Input tree file:          ", tree_file)
message("## Output figure file:       ", figure_file, "\n")



# PLOT THE TREE ----------------------------------------------------------------
## Read the tree file and prep the tree
tree <- read.tree(tree_file)
nseqs <- length(tree$tip.label)

## Make the plot
p <- ggtree(tree, layout = "rectangular") +
        geom_tiplab(size = 2.5) +
        geom_treescale(x = 0, y = nseqs - 1, color = "grey50") +
        geom_rootedge(rootedge = 0.005) +
        theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm")) +
        coord_cartesian(clip = "off")

ggsave(figure_file, p)

# WRAP UP ----------------------------------------------------------------------
message("\n## Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n## Done with script ggtree.R")
Sys.time()