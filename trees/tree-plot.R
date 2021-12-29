#!/usr/bin/env Rscript

# SBATCH --account=PAS0471
# SBATCH --time=5
# SBATCH --output=slurm-tree-plot-%j.out

# SET-UP -----------------------------------------------------------------------
## Command-line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1] # Input tree file
alignment_file <- args[2] # Input MSA file
annot_file <- args[3] # Input processed BLAST output file
fig_file <- args[4] # Output file with plot
if (length(args) >= 5) show_strain <- as.logical(args[5]) # TRUE/FALSE - whether to include strain info in tree tip labels
if (length(args) >= 6) msa_offset <- as.numeric(args[6]) # Offset of MSA in plot; use 0 for automatic determination
if (length(args) > 6) stop("Error: You provided more than 6 arguments")

if (!exists("show_strain")) show_strain <- FALSE
if (!exists("msa_offset")) msa_offset <- "auto"
if (is.na(msa_offset) | msa_offset == 0) msa_offset <- "auto"

# tree_file <- "results/tree/vir_ass/aedes/zvirus_mafft_blast.tre"
# alignment_file <- "results/tree/vir_ass/aedes/zvirus_mafft_blasthits_aligned.aln"
# annot_file <- "results/blast/vir_ass/aedes/AeJap_blast.proc"
# fig_file <- "results/tree/vir_ass/aedes/zvirus_tree.png"
# show_strain <- FALSE
# msa_offset <- 0

## Report
message("## Starting script tree-plot.R")
message(Sys.time())
message("## Tree file: ", tree_file)
message("## Alignment file: ", alignment_file)
message("## Annotation file: ", annot_file)
message("## Tree figure file: ", fig_file)
message("## Show strain info: ", show_strain)
message("## MSA offset in plot: ", msa_offset)
message("--------------------------------\n\n")

## Test
stopifnot(file.exists(tree_file))
stopifnot(file.exists(alignment_file))
stopifnot(file.exists(annot_file))

## Load (and install, if necessary) packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse", "here", "ape", "ggtree")
pacman::p_load(char = packages, install = TRUE)

## Source script with functions
source(here("scripts/tree-plot_fun.R"))

## Process args
fig_dir <- dirname(fig_file)
alt_msa_dir <- file.path(fig_dir, "alt_msa")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(alt_msa_dir)) dir.create(alt_msa_dir, recursive = TRUE)


# CREATE TREE ------------------------------------------------------------------
## Read tree
tree <- read.tree(tree_file)

## Read and parse annotation file
annot <- prep_annot(annot_file, tree$tip.label, show_strain)

## Plot the tree
plot_tree(tree, annot, alignment_file, fig_file, msa_offset)

## Plot the tree w/ several different settings for the MSA offset
msa_offsets <- seq(0.05, 0.6, by = 0.05)
for (msa_offset in msa_offsets) {
    fig_file_msa_base <- sub(
        ".png", paste0("_msa", msa_offset, ".png"),
        basename(fig_file)
    )
    fig_file_msa <- file.path(alt_msa_dir, fig_file_msa_base)
    plot_tree(tree, annot, alignment_file, fig_file_msa, msa_offset)
}

# WRAP UP ----------------------------------------------------------------------
if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

message("## Listing output file:")
system(paste("ls -lh", fig_file))
message("## Done with script tree-plot.R")
message(Sys.time())