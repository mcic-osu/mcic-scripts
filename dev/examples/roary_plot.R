#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --job-name=roary-plot
#SBATCH --output=slurm-roary_plot-%j.out

# SET-UP -----------------------------------------------------------------------
packages <- c("pagoo", "argparse")
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(BiocManager)) install.packages("pagoo", repos = rep, lib = lib)
} )
pacman::p_load(char = packages, install = TRUE)

# Parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--infile",
                    type = "character", required = TRUE,
                    help = "Roary's 'gene_presence_absence.csv' output file (REQUIRED)")
parser$add_argument("-o", "--figure",
                    type = "character", required = FALSE, default = NULL,
                    help = "Output figure file in '.png' format
                            [default: 'pangenome_rarefaction.png' in same dir as input file]")
args <- parser$parse_args()
infile <- args$infile
figure_file <- args$figure

# Define the output file name, if needed
if (is.null(figure_file)) {
  figure_file <- file.path(dirname(infile), "pangenome_rarefaction.png")
}

# Create outdir if needed
outdir <- dirname(figure_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Report
message("\n# Starting script roary_plot.R")
message("# Input file:                         ", infile)
message("# Output figure file:                 ", figure_file)
message()

# PLOT -------------------------------------------------------------------------
# Load the input
pg <- roary_2_pagoo(roary_outfile)

# Create the pangenome rarefaction graph (https://iferres.github.io/pagoo/articles/Methods_Plots.html)
pg$gg_curves() + 
  geom_point() + 
  facet_wrap(~Category, scales = 'free_y') +
  scale_color_brewer(palette = "Accent", guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave(figure_file, width = 7.5, height = 6, dpi = "retina")

# WRAP UP ----------------------------------------------------------------------
message("\n# Listing the output file:")
system(paste("ls -lh", figure_file))
message("\n# Done with script roary_plot.R")
Sys.time()
