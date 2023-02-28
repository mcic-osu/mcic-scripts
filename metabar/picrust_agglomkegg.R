#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=1
#SBATCH --time=15
#SBATCH --output=slurm-agglom-pathway-%j.out

# SET-UP -----------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse")
pacman::p_load(char = packages)

# Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ko_mat_in <- args[1]
kegg_map_in <- args[2]
mod_mat_out <- args[3]
path_mat_out <- args[4]

# Parameters for interactive testing
# ko_mat_in <- here("results/picrust/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
# kegg_map_in <- here("results/picrust/kegg_map.txt")
# mod_mat_out <- here("results/picrust/agglom_pathway/kegg_module_counts.txt")
# path_mat_out <- here("results/picrust/agglom_pathway/kegg_pathway_counts.txt")

# Create output dir if needed
outdir_mod <- dirname(mod_mat_out)
outdir_path <- dirname(path_mat_out)
dir.create(outdir_mod, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_path, recursive = TRUE, showWarnings = FALSE)

# Report
message("\n# Starting script agglom_pathway.R")
Sys.time()
message("# KEGG count matrix (input):           ", ko_mat_in)
message("# KEGG map/lookup (input) :            ", kegg_map_in)
message("# Module count matrix (output):        ", mod_mat_out)
message("# Pathway count matrix (output):       ", path_mat_out)


# PREPARE THE INPUT FILES ------------------------------------------------------
# KEGG lookup table
kegg_map <- read_tsv(kegg_map_in, show_col_types = FALSE) %>%
  dplyr::select(-symbol)

# KO count matrix
ko_mat <- read_tsv(ko_mat_in, show_col_types = FALSE) %>% as.data.frame()
row.names(ko_mat) <- ko_mat[, 1]     # First column => row names
ko_mat <- as.matrix(ko_mat[, -1])
ko_mat <- round(ko_mat)              # Round to integers or aldex will complain
ko_mat <- ko_mat[, !grepl("^WW", colnames(ko_mat))] # Remove wastewater samples


# MATRIX BY KEGG MODULE --------------------------------------------------------
# Get all modules
modules <- kegg_map %>%
  filter(KO_id %in% rownames(ko_mat)) %>%
  pull(pw_mod) %>%
  .[grepl("^M", .)] %>%
  unique()

# Create abundance matrix w/ one row per pathway, summing KOs in each module
mod_mat <- do.call(rbind,
                   lapply(X = modules, FUN = getmod, ko_mat, meta, kegg_map)) %>%
  as.data.frame() %>% 
  rownames_to_column("kegg_module")

write_tsv(mod_mat, mod_mat_out)


# MATRIX BY KEGG PATHWAY -------------------------------------------------------
# Get all pathways
pathways <- kegg_map %>%
  filter(KO_id %in% rownames(ko_mat)) %>%
  pull(pw_mod) %>%
  .[grepl("^map", .)] %>%
  unique()

# Create abundance matrix w/ one row per pathway, summing KOs in each pathway
path_mat <- do.call(rbind,
                    lapply(X = pathways, FUN = getmod, ko_mat, meta, kegg_map)) %>%
  as.data.frame() %>% 
  rownames_to_column("kegg_module")

write_tsv(path_mat, path_mat_out)


# WRAP UP ----------------------------------------------------------------------
message("\n# Listing the output files:\n")
system(paste("ls -lh", mod_mat_out))
system(paste("ls -lh", path_mat_out))

message("\n# Done with script agglom_pathway.R")
Sys.time()
sessionInfo()
