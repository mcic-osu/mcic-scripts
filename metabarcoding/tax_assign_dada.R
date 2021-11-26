#!/usr/bin/env Rscript

# SET-UP -----------------------------------------------------------------------
## Load packages
if(!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("BiocManager", "dada2", "DECIPHER", "tidyverse")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
seqtab_rds <- args[1]
taxa_rds <- args[2]
n_cores <- as.integer(args[3])

# seqtab_rds <- "results/ASV/main/seqtab.rds"
# taxa_rds <- "results/taxonomy/taxa_dada.rds"
# n_cores <- 8

## Create output dir if needed
outdir <- dirname(taxa_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Other variables/constants
tax_URL <- "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
tax_file <- file.path(outdir, basename(tax_URL))
species_URL <- "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"
species_file <- file.path(outdir, basename(species_URL))

qc_file <- file.path(outdir, "tax_prop_assigned_dada.txt")
plot_file <- file.path(outdir, "tax_prop_assigned_dada.png")

tax_levels <- c("domain", "phylum", "class", "order", "family", "genus", "species")

## Report
cat("\n## Starting script tax_assign_dada.R\n")
Sys.time()
cat("## Sequence table RDS file (input):", seqtab_rds, "\n")
cat("## Taxa RDS file (output):", taxa_rds, "\n")
cat("## Proportion-assigned QC file (output):", qc_file, "\n")
cat("## Number of cores:", n_cores, "\n\n")
cat("## Taxonomic assignment file (downloaded input):", tax_file, "\n")
cat("## Species assignment file (downloaded input):", species_file, "\n\n")

## Function to get the proportion of ASVs assigned to taxa
qc_tax <- function(taxa, tax_levels) {
    prop <- apply(taxa, 2,
                  function(x) round(length(which(!is.na(x))) / nrow(taxa), 4))
    
    prop <- data.frame(prop) %>%
        rownames_to_column("tax_level") %>%
        mutate(tax_level = factor(tax_level, levels = tax_levels))

    return(prop)
}


# PREPARE INPUT DATA -----------------------------------------------------------
## Create a DNAStringSet from the ASVs
seqtab <- readRDS(seqtab_rds)
dna <- DNAStringSet(getSequences(seqtab))


# DADA2 TAX. ASSIGNMENT --------------------------------------------------------
## Get and load DADA training set
## (Check for an up-to-date version at <https://benjjneb.github.io/dada2/training.html>)
if (!file.exists(tax_file)) download.file(url = tax_URL, destfile = tax_file)
if (!file.exists(species_file)) download.file(url = species_URL, destfile = species_file)

## Assign taxonomy
cat("\n## Now assigning taxonomy...\n")
taxa <- assignTaxonomy(seqtab, tax_file, multithread = n_cores)

cat("\n## Now adding species-level assignments...\n")
taxa <- addSpecies(taxa, species_file)
colnames(taxa) <- tax_levels

## Save RDS file
saveRDS(taxa, taxa_rds)


# QC TAX. ASSIGNMENTS ----------------------------------------------------------
## Create df with proportion assigned
prop_assigned <- qc_tax(taxa, tax_levels = tax_levels)
write.table(prop_assigned, qc_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n## Proportion of ASVs assigned to different taxonomic levels - DADA:\n")
print(prop_assigned)

## Create barplot
p <- ggplot(prop_assigned) +
    geom_col(aes(x = tax_level, y = prop, fill = tax_level),
             color = "grey20") +
    scale_fill_brewer(palette = "Greens") +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(y = "Proportion of ASVs assigned", x = NULL) +
    guides(fill = "none") +
    theme_bw(base_size = 14)
ggsave(plot_file, p, width = 7, height = 7)


# WRAP UP ----------------------------------------------------------------------
## Report
cat("\n## Listing output files:\n")
system(paste("ls -lh", taxa_rds))
system(paste("ls -lh", qc_file))
system(paste("ls -lh", plot_file))

cat("\n## Done with script tax_assign_dada.R\n")
Sys.time()
