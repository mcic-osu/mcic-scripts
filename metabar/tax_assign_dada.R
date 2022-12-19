#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-tax_assign-%j.out


# SET-UP -----------------------------------------------------------------------
# Packages
packages <- c("BiocManager", "dada2", "DECIPHER", "tidyverse")

# Constants
TAX_LEVELS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Defaults
ref_url <- "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
species_url <- "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"

# Set nr of threads
n_threads <- as.integer(system("echo $SLURM_CPUS_PER_TASK", intern = TRUE))

# Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)
parser <- ArgumentParser()
parser$add_argument("-i", "--seq_file",
                    type = "character", required = TRUE,
                    help = "Input file (sequence table RDS or Qiime2 QZA) (REQUIRED)")
parser$add_argument("-o", "--taxa_out",
                    type = "character", required = TRUE,
                    help = "Output file (taxa RDS file) (REQUIRED)")
parser$add_argument("-r", "--ref_file",
                    type = "character",
                    default = NULL,
                    help = "Taxonomic reference file [default: download Silva v138.1]")
parser$add_argument("--species_file",
                    type = "character",
                    default = NULL,
                    help = "Taxonomic reference file [default download Silva v138.1 species file]")
parser$add_argument("--add_species",
                    type = "logical",
                    default = TRUE,
                    help = "Add species-level taxonomy")
args <- parser$parse_args()

seq_file <- args$seq_file
taxa_rds <- args$taxa_out
ref_file <- args$ref_file
species_file <- args$species_file
add_species <- args$add_species

# Define output files
outdir <- dirname(taxa_rds)
outdir_ref <- file.path(outdir, "tax_ref")
outdir_qc <- file.path(outdir, "qc")
qc_file <- file.path(outdir_qc, "tax_prop_assigned_dada.tsv")
plot_file <- file.path(outdir_qc, "tax_prop_assigned_dada.png")

# Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(char = packages, repos = "https://cran.rstudio.com/")

# Create output dir if needed
dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_ref, recursive = TRUE, showWarnings = FALSE)

# Report
message()
message("# ====================================================================")
message("#               STARTING SCRIPT TAX_ASSIGN_DADA.R")
message("# ====================================================================")
Sys.time()
message()
message("# Input file with sequences:                     ", seq_file)
message("# Taxa RDS file (output):                        ", taxa_rds)
message("# Proportion-assigned QC file (output):          ", qc_file)
message("# Number of cores:                               ", n_threads)
message("# Taxonomic assignment file (downloaded input):  ", ref_file)
message("# Species assignment file (downloaded input):    ", species_file)
message("# ====================================================================")
message()


# FUNCTIONS --------------------------------------------------------------------
# Function to get the proportion of ASVs assigned to taxa
qc_tax <- function(taxa_df, tax_levels = NULL) {
    
    n <- apply(taxa_df, 2, function(x) length(which(!is.na(x))))
    prop <- round(n / nrow(taxa_df), 4)
    
    if (is.null(tax_levels)) tax_levels <- names(prop)
    colnames(taxa) <- tax_levels
    
    data.frame(n, prop) %>%
        rownames_to_column("tax_level") %>%
        mutate(tax_level = factor(tax_level, levels = tax_levels))
}


# PREPARE REFERENCE SEQUENCES --------------------------------------------------
if (is.null(ref_file)) {
    ref_file <- file.path(outdir_ref, basename(ref_url))
    if (!file.exists(ref_file)) download.file(url = ref_url, destfile = ref_file)
}
if (is.null(species_file)) {
    species_file <- file.path(outdir_ref, basename(species_url))
    if (!file.exists(species_file)) download.file(url = species_url, destfile = species_file)
}

# PREPARE INPUT SEQUENCES ------------------------------------------------------
if (grepl("\\.rds", seq_file, ignore.case = TRUE)) {
    seqs <- readRDS(seq_file)
    #dna <- DNAStringSet(getSequences(seq_file))
} else if (grepl("\\.qzv", seq_file, ignore.case = TRUE)) {
    seqs <- read_qza(seq_artif)$data
}


# DADA2 TAX. ASSIGNMENT --------------------------------------------------------
message("# Now assigning taxonomy...")
taxa <- assignTaxonomy(seqs, ref_file, multithread = n_threads)

if (add_species == TRUE) {
    message("\n# Adding species-level assignments...")
    taxa <- addSpecies(taxa, species_file)
}

# Save RDS file
saveRDS(taxa, taxa_rds)


# QC TAX. ASSIGNMENTS ----------------------------------------------------------
# Create df with proportion assigned
prop_assigned <- qc_tax(taxa_df = taxa, tax_levels = TAX_LEVELS)
write_tsv(prop_assigned, qc_file)

message("\n# Prop. ASVs assigned to taxonomy:")
print(prop_assigned)

# Create barplot
p <- ggplot(prop_assigned) +
    geom_col(aes(x = tax_level, y = prop, fill = tax_level), color = "grey20") +
    scale_fill_brewer(palette = "Greens", direction = -1) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(y = "Proportion of ASVs assigned", x = NULL) +
    guides(fill = "none") +
    theme_bw(base_size = 14)
ggsave(plot_file, p, width = 7, height = 7)


# WRAP UP ----------------------------------------------------------------------
message("\n# Listing output files:")
system(paste("ls -lh", taxa_rds))
system(paste("ls -lh", qc_file))
system(paste("ls -lh", plot_file))

message("\n# Done with script tax_assign_dada.R")
Sys.time()
message()
