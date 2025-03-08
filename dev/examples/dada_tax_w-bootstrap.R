#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --job-name=dada-assign-tax
#SBATCH --output=slurm-dada-assign-tax-%j.out

# Load packages
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))

# Parse command-line args
args <- commandArgs(trailingOnly = TRUE)
ref_file <- args[1]  # Reference FASTA file, formatted for dada assignTaxonomy()
asv_file <- args[2]  # ASV FASTA file
outdir <- args[3]

# Report
message("\n# Starting script dada_eukaryome_assign-tax.R")
Sys.time()
message("# Reference FASTA file:  ", ref_file)
message("# Reference FASTA file:  ", asv_file)
message("# Output dir:            ", outdir, "\n")

# Define the output files
outdir_log <- file.path(outdir, "logs")
dir.create(outdir_log, showWarnings = FALSE, recursive = TRUE)
tax_df_file <- file.path(outdir, "eukaryome_dada_tax.tsv")

# Read the ASV file
asv_seqs <- getSequences(asv_file)

# Assign taxonomy
message("# Now assigning taxonomy...")
taxonomy <- assignTaxonomy(
  seqs = asv_seqs,
  refFasta = ref_file,
  tryRC = TRUE,
  minBoot = 50,
  outputBootstraps = TRUE,
  multithread = 8,
  )

message("# Now processing the output...")
# ASV seqs dataframe
asv_df <- tibble(seq = asv_seqs, ASV = names(asv_seqs))

write_tsv(data.frame(taxonomy$tax), file.path(outdir, "tax.tsv"))
write_tsv(data.frame(taxonomy$boot), file.path(outdir, "boot.tsv"))

# Reformat the taxonomy output
tax_df <- data.frame(taxonomy$tax) |>
  rownames_to_column("seq") |>
  full_join(asv_df, by = "seq") |>
  select(-seq) |> 
  tibble() |>
  relocate(ASV)

# Process the bootstrap support
boot_df <- data.frame(taxonomy$boot) |>
  rownames_to_column("seq") |>
  full_join(asv_df, by = "seq") |>
  select(-seq) |> 
  tibble() |>
  relocate(ASV)

# Join the dataframes
tax_boot_df <- full_join(tax_df, boot_df, by = "ASV", suffix = c("", "_boot"))

# Write the output files
write_tsv(tax_boot_df, tax_df_file)

# Report
message("\n# Done with script dada_eukaryome_assign.R")
message("# Listing the output file:")
system(paste("ls -lh", tax_df_file))
Sys.time()
