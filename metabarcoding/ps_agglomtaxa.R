#!/usr/bin/env Rscript

#SBATCH --account=PAS0471 # nolint
#SBATCH --output=slurm-ps_split-%j.out # nolint

# SETUP ------------------------------------------------------------------------
## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c(
    "BiocManager", "tidyverse", "phyloseq", "decontam",
    "microbiome", "biomformat"
)
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ps_in <- args[1]
outdir <- args[2]

## Define output files
file_id <- sub(".rds", "", basename(ps_in))

outfile_phylum <- file.path(outdir, paste0(file_id, "_phylum.rds"))
outfile_class <- file.path(outdir, paste0(file_id, "_class.rds"))
outfile_order <- file.path(outdir, paste0(file_id, "_order.rds"))
outfile_family <- file.path(outdir, paste0(file_id, "_family.rds"))
outfile_genus <- file.path(outdir, paste0(file_id, "_genus.rds"))
outfile_species <- file.path(outdir, paste0(file_id, "_species.rds"))

## Report
cat("\n## Starting script ps_split.R\n")
Sys.time()
cat("## Input phyloseq RDS file:             ", ps_in, "\n")
cat("## Output dir:                          ", outdir, "\n")
cat("-----------------------------\n\n")


# CREATE PHYLOSEQ OBJECTS AGGLOMERATED BY TAXRANK ------------------------------
## Read input file
ps <- readRDS(ps_in)

## Agglomerate
ps_phylum <- tax_glom(ps, taxrank = "phylum")
ps_class <- tax_glom(ps, taxrank = "class")
ps_order <- tax_glom(ps, taxrank = "order")
ps_family <- tax_glom(ps, taxrank = "family")
ps_genus <- tax_glom(ps, taxrank = "genus")
ps_species <- tax_glom(ps, taxrank = "species")

saveRDS(ps_phylum, outfile_phylum)
saveRDS(ps_class, outfile_class)
saveRDS(ps_order, outfile_order)
saveRDS(ps_family, outfile_family)
saveRDS(ps_genus, outfile_genus)
saveRDS(ps_species, outfile_species)


# WRAP UP ----------------------------------------------------------------------
cat("\n\n## Listing output files:\n")
system(paste("ls -lh", outdir))

cat("\n## Done with script ps_split.R\n")
Sys.time()