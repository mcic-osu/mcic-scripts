#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=10
#SBATCH --output=slurm-tree-plot-%j.out

# SET-UP -----------------------------------------------------------------------
# Load (and install, if necessary) packages
#chooseCRANmirror(ind = 1)
if (!"pacman" %in% installed.packages()) {
    install.packages("pacman",
                    repos = "https://cloud.r-project.org/",
                    lib = Sys.getenv("R_LIBS_USER"))
}
packages <- c("tidyverse", "here", "ape", "BiocManager", "ggtree", "argparse")
pacman::p_load(char = packages, install = TRUE,
                repos = "https://cloud.r-project.org/",
                lib = Sys.getenv("R_LIBS_USER"))

# Other scripts
fun_script <- "mcic-scripts/trees/tree-plot_fun.R"

# Parse options
parser <- ArgumentParser() # create parser object
parser$add_argument("-i", "--tree",
                    help = "Input tree file (REQUIRED)")
parser$add_argument("-o", "--figure",
                    help = "Output figure (REQUIRED)")
parser$add_argument("-a", "--aln",
                    help = "Input alignment file (required if --show_msa=TRUE)")
parser$add_argument("-n", "--annot", default = NULL,
                    help = "Input annotation file")
parser$add_argument("-l", "--label_col1", type = "character", default = NULL,
                    help = "Name of 1st annotation file column to show as tip label")
parser$add_argument("-L", "--label_col2", type = "character", default = NULL,
                    help = "Name of 2nd annotation file column to show as tip label")
parser$add_argument("--label_col3", type = "character", default = NULL,
                    help = "Name of 3rd annotation file column to show as tip label")
parser$add_argument("-c", "--color_col", type = "character", default = NULL,
                    help = "Name of annotation file column to color tip labels by")
parser$add_argument("--blast", action = "store_true", default = FALSE,
                    help = "Annotation file is BLAST output")
parser$add_argument("--plot_msa", type = "logical", default = FALSE,
                    help = "Show MSA")
parser$add_argument("--msa_offset", default = "auto",
                    help = "Offset for MSA")
parser$add_argument("--show_strain", type = "logical", default = FALSE,
                    help = "Show strain")
args <- parser$parse_args()

if (! args$msa_offset %in% c("auto", "textlen")) {
    args$msa_offset <- as.numeric(args$msa_offset)
}

# Test input
stopifnot(file.exists(args$tree))
if (args$plot_msa == TRUE) stopifnot(file.exists(args$aln))
if (!is.null(args$annot)) stopifnot(file.exists(args$annot))

# Source script with functions
source(here(fun_script))

# Process args
fig_dir <- dirname(args$fig)
alt_msa_dir <- file.path(fig_dir, "alt_msa")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(alt_msa_dir)) dir.create(alt_msa_dir, recursive = TRUE)

if (!is.null(args$label_col1)) {
    if (is.null(args$label_col3)) {
        if (is.null(args$label_col2)) {
            label_cols <- args$label_col1
        } else {
            label_cols <- c(args$label_col1, args$label_col2)
        }
    } else {
        label_cols <- c(args$label_col1, args$label_col2, args$label_col3)
    }
} else {
    label_cols <- NULL
}

# Report
message("\n# Starting script tree-plot.R")
message(Sys.time())
message()
message("# Input tree file:            ", args$tree)
if (!is.null(args$annot)) message("# Annotation file:            ", args$annot)
message("# Plot MSA besides the tree:  ", args$plot_msa)
if (args$plot_msa == TRUE) message("# Alignment file:              ", args$aln)
message("# Output figure:              ", args$figure)
if (!is.null(label_cols)) message("# Label columns:              ", label_cols)
if (args$show_strain == TRUE) message("# Show strain info:            ", args$show_strain)
if (args$plot_msa == TRUE) message("# MSA offset in plot:          ", args$msa_offset)
message("----------------------\n")

# CREATE TREE ------------------------------------------------------------------
# Read the tree file
tree <- read.tree(args$tree)

# Read and parse annotation file
if (!is.null(args$annot)) {
    if (args$blast == TRUE) {
        annot <- prep_annot_blast(args$annot, tree$tip.label, show_strain)
    } else {
        annot <- read_tsv(args$annot, show_col_types = FALSE)
    }
} else {
    annot <- NULL
}

# Plot the tree
plot_tree(tree = tree,
          alignment = args$aln,
          fig_file = args$fig,
          annot = annot,
          label_columns = label_cols,
          color_column = args$color_col,
          plot_msa = args$plot_msa,
          msa_offset = args$msa_offset)

# WRAP UP ----------------------------------------------------------------------
if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

message("# Listing output file:")
system(paste("ls -lh", args$fig))
message()
message("# Done with script tree-plot.R")
message(Sys.time())
message()
