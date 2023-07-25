#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --job-name=ggtimetree
#SBATCH --output=slurm-ggtimetree-%j.out

#? From an input tree file, this script will plot a time-tree with ggtree
#? (All tree file formats should be supported)

#? Load the Conda environment as follows to run this script without it needing to install R packages:
#? module load miniconda3 && source activate /fs/ess/PAS0471/jelmer/conda/r_tree

# SET-UP -----------------------------------------------------------------------
# Load packages
packages <- c("argparse", "ggtree", "ggplot2")
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
parser$add_argument("--dates",
                    type = "character", required = FALSE, default = NULL,
                    help = "File with dates (use either --dates or --mrsd)
                            Dates should be in the second column")
parser$add_argument("--header",
                    action = "store_true", required = FALSE, default = FALSE,
                    help = "Dates file has a header [default: false]")
parser$add_argument("--mrsd",
                    type = "character", required = FALSE, default = NULL,
                    help = "Most recent sampling date in YYYY-MM-DD format (use either --dates or --mrsd)")
args <- parser$parse_args()

tree_file <- args$tree
figure_file <- args$figure
dates_file <- args$dates
has_header <- args$header
mrsd <- args$mrsd

# Define the output file name, if needed
if (is.null(figure_file)) {
  figure_file <- file.path(
    dirname(tree_file),
    paste0(tools::file_path_sans_ext(basename(tree_file)), ".png")
    )
}

# Get the most recent sampling date from the dates file
if (!is.null(dates_file)) {
  dates <- read.delim(dates_file, header = has_header)
  mrsd <- sort(dates[[2]])[length((dates[[2]]))]
}
if (is.null(mrsd)) stop("No most recent sample date, use --dates or --mrsd")

# Report
message("\n# Starting script ggtimetree.R")
message("# Input tree file:          ", tree_file)
if (!is.null(dates_file)) message("# Input dates file:         ", dates_file)
message("# Most recent sample date:  ", mrsd)
message("# Output figure file:       ", figure_file, "\n")


# PLOT THE TREE ----------------------------------------------------------------
# Read the tree
tree_file_ext <- tools::file_ext(tree_file)
if (tree_file_ext == "nex") {
  tree <- treeio::read.nexus(tree_file)
} else {
  tree <- treeio::read.tree(tree_file)
}

# Make the plot
p <- ggtree(tree, layout = "rectangular", mrsd = mrsd) +
  geom_tiplab(size = 4, color = "grey50") +
  geom_rootedge(rootedge = sum(tree$edge.length) / 50) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(n.breaks = 8) +
  theme_tree2() +
  theme(plot.margin = margin(0.2, 3, 0.2, 0.2, "cm"),
        axis.text.x = element_text(size = 13, color = "grey50"),
        axis.line.x.bottom = element_line(color = "grey50"),
        axis.ticks.x.bottom = element_blank(),
        panel.grid.major.x = element_line(linetype = "longdash"))

# Save the plot to file
ggsave(figure_file, p, width = 6, height = 6, dpi = "retina")


# WRAP UP ----------------------------------------------------------------------
message("# Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n# Done with script ggtimetree.R")
Sys.time()
