#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-kegg_map-%j.out
#SBATCH --time=12:00:00

# Script to create a KEGG 'map'/lookup
# Input is a file that contains KO numbers in the first column
# The script will look up KEGG modules and pathways that each KO number belongs to

# Example command to run the script:
#> ko_list=/fs/ess/PAS0471/jelmer/refdata/kegg/ko_list
#> ko_map=/fs/ess/PAS0471/jelmer/refdata/kegg/ko_map.tsv
#> ml R/4.2.1-gnu11.2 && sbatch mcic-scripts/annot/kegg_map.R "$ko_list" "$ko_map"


# SET-UP -----------------------------------------------------------------------
# Install/load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c(
    "tidyverse",
    "here",
    "KEGGREST",
    "clusterProfiler"
)
pacman::p_load(char = packages)

# Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ko_file <- args[1]             # Input file with KO numbers
kegg_map_file <- args[2]       # Output: KO-to-pathway+module mappings

# Example arguments
# ko_file <- "/fs/ess/PAS0471/jelmer/refdata/kegg/ko_list"
# kegg_map_file <- "/fs/ess/PAS0471/jelmer/refdata/kegg/ko_map.tsv"

# Create output dir if needed
outdir <- dirname(kegg_map_file)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Report
message("# Starting script kegg_map.R")
Sys.time()
message("# Input KO file:               ", ko_file)
message("# Output KEGG map:             ", kegg_map_file)
message("==============================================\n")

# Check input
stopifnot(file.exists(ko_file))


# FUNCTIONS --------------------------------------------------------------------
# Function to get a KEGG pathway (ko-ID) associated with a KEGG K-term
ko_info <- function(ko_id, outdir) {
  cat("# KO ID:", ko_id, " ")
  tryCatch(
    {
      keggres <- keggGet(ko_id)

      if (!is.null(keggres[[1]]$PATHWAY)) pw <- names(keggres[[1]]$PATHWAY) else pw <- NA
      if (!is.null(keggres[[1]]$PATHWAY)) pw_descr <- keggres[[1]]$PATHWAY else pw_descr <- NA
      if (!is.null(keggres[[1]]$MODULE)) mod <- names(keggres[[1]]$MODULE) else mod <- NA
      if (!is.null(keggres[[1]]$MODULE)) mod_descr <- keggres[[1]]$MODULE else mod_descr <- NA

      pathway_df <- data.frame(
        KO_id = ko_id,
        symbol = ifelse(is.null(keggres[[1]]$SYMBOL), NA, keggres[[1]]$SYMBOL),
        name = ifelse(is.null(keggres[[1]]$NAME), NA, keggres[[1]]$NAME),
        pw_mod = pw,
        pw_mod_descr = pw_descr,
        row.names = NULL
      )
      cat("// n pathways:", nrow(pathway_df))

      if (!is.na(mod[1])) {
        module_df <- data.frame(
            KO_id = ko_id,
            symbol = ifelse(is.null(keggres[[1]]$SYMBOL), NA, keggres[[1]]$SYMBOL),
            name = ifelse(is.null(keggres[[1]]$NAME), NA, keggres[[1]]$NAME),
            pw_mod = mod,
            pw_mod_descr = mod_descr,
            row.names = NULL
        )
        kegg_df <- bind_rows(pathway_df, module_df)
        cat(" // n modules:", nrow(module_df), "\n")
      } else {
        kegg_df <- pathway_df
        cat("\n")
      }
      if (nrow(kegg_df) > 0) return(kegg_df) else return(NULL)
    },
    error = function(cond) {
      message("keggGet failure for ", cond)
      return(NULL)
    }
  )
}

# CREATE KEGG MAP --------------------------------------------------------------
# Read input files
ko_terms <- read_tsv(ko_file, show_col_types = FALSE) %>% pull(1)

# Get KEGG pathways associated with the different K-terms
message("\n# Now creating the KEGG map...")
kegg_map <- map_dfr(.x = ko_terms, .f = ko_info, outdir_terms)

# Write files
write_tsv(kegg_map, kegg_map_file)


# WRAP UP ----------------------------------------------------------------------
message("\n==============================================")
message("# Listing the output file:")
system(paste("ls -lh", kegg_map_file))

message("\n# Done with script kegg_map.R")
Sys.time()
sessionInfo()
