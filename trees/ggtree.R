#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --job-name=ggtree
#SBATCH --output=slurm-ggtree-%j.out

#TODO - Add node labels (bootstrap)
#p <- p +
#  geom_text2(aes(subset = !isTip, label = label),
#             color = "grey50", nudge_y = 0.4, nudge_x = -0.05, size = 3)

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
parser$add_argument("--root",
                    type = "character", default = NULL,
                    help = "ID of sample that should be the root of the tree: tree will be rerooted")
parser$add_argument("--annot",
                    type = "character", default = NULL,
                    help = "Input annotation/metadata file")
parser$add_argument("--layout",
                    type = "character", default = "rectangular",
                    help = "Tree layout")
parser$add_argument("--color_column",
                    type = "character", default = NULL,
                    help = "Name of annotation file column to color tip labels by")
parser$add_argument("--tiplab_column",
                    type = "character", default = NULL,
                    help = "Name of annotation file column to use as tip labels (instead of labels in the tree file)")
parser$add_argument("--right_margin",
                    type = "numeric", default = 2,
                    help = "Size of the plot's right margin: extend to avoid truncated tip labels, etc.")

args <- parser$parse_args()
tree_file <- args$tree
figure_file <- args$figure
annot_file <- args$annot
layout <- args$layout
color_column <- args$color_column
tiplab_column <- args$tiplab_column
root <- args$root
right_margin <- args$right_margin

# Define the output file name, if needed
if (is.null(figure_file)) {
  figure_file <- file.path(
    dirname(tree_file),
    paste0(tools::file_path_sans_ext(basename(tree_file)), ".png")
  )
}

# Report
message("\n# Starting script ggtree.R")
message("# Input tree file:                         ", tree_file)
message("# Output figure file:                      ", figure_file)
if (!is.null(annot_file)) message("# Annotation/metadata file:                ", annot_file)
if (!is.null(tiplab_column)) message("# Metadata column for tip labels:          ",      tiplab_column)
if (!is.null(color_column)) message("# Metadata column for colors:              ", color_column)
if (!is.null(root)) message("# ID of sample that should be the root:    ", root)
message()


# PLOT PREP --------------------------------------------------------------------
# Read the tree file and prep the tree
tree <- read.tree(tree_file)
nseqs <- length(tree$tip.label)

# Read the annotation
if (!is.null(annot_file)) {
  annot <- read.delim(annot_file)
  if (!is.null(tiplab_column)) {
    annot$tiplab <- annot[[tiplab_column]]
    tiplab_column <- "tiplab"
  }
  message("\n# Showing the first few lines of the annotation dataframe:")
  print(head(annot))
  cat("\n")

  # Check that IDs for annotation & tree match
  message("\n# Showing the tree's tip labels:")
  print(tree$tip.label)

  message("\n# Are all tip labels in the annotation df?")
  print(all(tree$tip.label %in% annot[[1]]))
}

# Tiplab size
if (nseqs < 50) LABEL_SIZE <- 2.5 else LABEL_SIZE <- 2

# Add root option
if (! is.null(root)) {
  message("\n# Rerooting the tree, using ", root, " as the root...")
  tree <- ape::root(tree, outgroup = root) 
}


# PLOT THE TREE ----------------------------------------------------------------
# Base tree
p <- ggtree(tree, layout = layout)

# Add annotation dataframe if provided
if (!is.null(annot_file)) p <- p %<+% annot

# Tip labels
if (!is.null(tiplab_column)) {
  p <- p + geom_tiplab(aes_string(color = color_column, label = tiplab_column),
                       align = TRUE, linesize = 0, size = LABEL_SIZE)
} else {
  p <- p + geom_tiplab(aes_string(color = color_column),
                       align = TRUE, linesize = 0, size = LABEL_SIZE)
}

# Make the plot
p <- p +
  geom_rootedge(rootedge = sum(tree$edge.length) / 50) +
  theme(plot.margin = margin(0.2, right_margin, 0.2, 0.2, "cm"),
        legend.box.spacing = unit(25, "pt"))
if(layout == "rectangular") p <- p + coord_cartesian(clip = "off")

# Save the plot to file
ggsave(figure_file, p, width = 6, height = 6, dpi = "retina")


# WRAP UP ----------------------------------------------------------------------
message("\n# Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n# Done with script ggtree.R")
Sys.time()
