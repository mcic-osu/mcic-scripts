#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-ps_agglomtaxa-%j.out


# SETUP ------------------------------------------------------------------------
# Packages
packages <- c("BiocManager", "phyloseq")

# Process command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--ps",
                    type = "character", required = TRUE,
                    help = "Input file (phyloseq object RDS) (REQUIRED)")
parser$add_argument("-o", "--outdir",
                    type = "character", required = TRUE,
                    help = "Output dir (REQUIRED)")
args <- parser$parse_args()
ps_in <- args$ps
outdir <- args$outdir

# Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(char = packages)

# Define output files
file_id <- sub(".rds", "", basename(ps_in))
outfile_phylum <- file.path(outdir, paste0(file_id, "_phylum.rds"))
outfile_class <- file.path(outdir, paste0(file_id, "_class.rds"))
outfile_order <- file.path(outdir, paste0(file_id, "_order.rds"))
outfile_family <- file.path(outdir, paste0(file_id, "_family.rds"))
outfile_genus <- file.path(outdir, paste0(file_id, "_genus.rds"))
outfile_species <- file.path(outdir, paste0(file_id, "_species.rds"))

# Report
message()
message("# ====================================================================")
message("#               STARTING SCRIPT PS_AGGLOMTAXA.R")
message("# ====================================================================")
Sys.time()
message("# Input phyloseq RDS file:             ", ps_in)
message("# Output dir:                          ", outdir)
message("# ====================================================================")
message()


# CREATE PHYLOSEQ OBJECTS AGGLOMERATED BY TAXRANK ------------------------------
# Read input file
ps <- readRDS(ps_in)

# Agglomerate
ps_phylum <- tax_glom(ps, taxrank = "phylum")
ps_class <- tax_glom(ps, taxrank = "class")
ps_order <- tax_glom(ps, taxrank = "order")
ps_family <- tax_glom(ps, taxrank = "family")
ps_genus <- tax_glom(ps, taxrank = "genus")
ps_species <- tax_glom(ps, taxrank = "species")

# Save output files
saveRDS(ps_phylum, outfile_phylum)
saveRDS(ps_class, outfile_class)
saveRDS(ps_order, outfile_order)
saveRDS(ps_family, outfile_family)
saveRDS(ps_genus, outfile_genus)
saveRDS(ps_species, outfile_species)


# WRAP UP ----------------------------------------------------------------------
message("\n# Listing output files:")
system(paste("ls -lh", outdir))

message("\n# Done with script ps_agglomtaxa.R")
Sys.time()
message()
sessionInfo()
message()
