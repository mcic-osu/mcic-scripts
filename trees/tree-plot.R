#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=5
#SBATCH --output=slurm-tree-plot-%j.out

# SET-UP -----------------------------------------------------------------------
## Load (and install, if necessary) packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse", "here", "ape", "ggtree", "argparse")
pacman::p_load(char = packages, install = TRUE)

## Other scripts
fun_script <- "mcic-scripts/trees/tree-plot_fun.R"

## Parse options
parser <- ArgumentParser() # create parser object
parser$add_argument("-t", "--tree",
                    help = "Input tree file (REQUIRED)")
parser$add_argument("-a", "--aln",
                    help = "Input alignment file (REQUIRED)")
parser$add_argument("-o", "--figure",
                    help = "Output figure (REQUIRED)")
parser$add_argument("-n", "--annot", default = NULL,
                    help = "Input annotation file")
parser$add_argument("-c", "--annot_col", type = "integer", default = 2,
                    help = "Annotation file column to show as tip label")
parser$add_argument("--blast", action = "store_true", default = FALSE,
                    help = "Annotation file is BLAST output")
parser$add_argument("--show_strain", type = "logical", default = FALSE,
                    help = "Show strain")
parser$add_argument("--msa_offset", default = "auto",
                    help = "Offset for MSA")
parser$add_argument("--try_msa_offsets", type = "logical", default = FALSE,
                    help = "Try many different MSA offsets")
args <- parser$parse_args()

if (! args$msa_offset %in% c("auto", "textlen")) {
    args$msa_offset <- as.numeric(args$msa_offset)
}

## Report
message("\n## Starting script tree-plot.R")
message(Sys.time())
message("## Tree file: ", args$tree)
message("## Alignment file: ", args$aln)
message("## Tree figure file: ", args$figure)
if (!is.null(args$annot)) message("## Annotation file: ", args$annot)
message()
message("## Show strain info: ", args$show_strain)
message("## MSA offset in plot: ", args$msa_offset)
message("## Try many MSA offsets: ", args$try_msa_offsets)
message("--------------------------------\n\n")

## Test
stopifnot(file.exists(args$tree))
stopifnot(file.exists(args$aln))
if (!is.null(args$annot)) stopifnot(file.exists(args$annot))

## Source script with functions
source(here(fun_script))

## Process args
fig_dir <- dirname(args$fig)
alt_msa_dir <- file.path(fig_dir, "alt_msa")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(alt_msa_dir)) dir.create(alt_msa_dir, recursive = TRUE)


# CREATE TREE ------------------------------------------------------------------
## Read the tree file
tree <- read.tree(args$tree)

## Read and parse annotation file
if (!is.null(args$annot)) {
    if (args$blast == TRUE) {
        annot <- prep_annot_blast(args$annot, tree$tip.label, show_strain)
    } else {
        annot <- prep_annot(args$annot, tree$tip.label)
    }
} else {
    annot <- NULL
}

## Plot the tree
plot_tree(tree, annot, args$aln, args$fig, args$msa_offset)

## Plot the tree w/ several different settings for the MSA offset
if (args$try_msa_offsets == TRUE) {
    for (offset in seq(0.05, 0.6, by = 0.05)) {
        fig <- here(alt_msa_dir,
                    sub(".png", paste0("_msa", offset, ".png"), basename(args$fig)))
        plot_tree(tree, annot, args$aln, fig, msa_offset)
    }
}

# WRAP UP ----------------------------------------------------------------------
if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

message("## Listing output file:")
system(paste("ls -lh", args$fig))
message("## Done with script tree-plot.R")
message(Sys.time())
message("\n")