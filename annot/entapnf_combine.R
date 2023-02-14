#!/usr/bin/env Rscript
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=entapnf_combine
#SBATCH --output=slurm-entapnf_combine-%j.out

# SET-UP -----------------------------------------------------------------------
# Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(tidyverse)) install.packages("tidyverse", repos = rep, lib = lib)
  if (!require(janitor)) install.packages("janitor", repos = rep, lib = lib)
} )
packages <- c("argparse", "tidyverse")
pacman::p_load(char = packages, install = TRUE)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--annot1",
                    type = "character", required = TRUE, default = NULL,
                    help = "Annotation file 1")
parser$add_argument("-I", "--annot2",
                    type = "character", required = TRUE, default = NULL,
                    help = "Annotation file 2")
parser$add_argument("-o", "--annot_out",
                    type = "character", required = TRUE, default = NULL,
                    help = "Output annotation file")
args <- parser$parse_args()

annot1_file <- args$annot1
annot2_file <- args$annot2
annot_out_file <- args$annot_out

#annot1_file <- "results/entap_proc_all/annotation.tsv"
#annot2_file <- "results/entap_proc_all_nr/annotation.tsv"

# Get and make the output dir
outdir <- dirname(annot_out_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Filenames for KEGG and GO output
go_file <- file.path(outdir, "go_map.tsv")
kegg_file <- file.path(outdir, "kegg_map.tsv")

# Report
message("\n# Starting script entapnf_combine.R")
Sys.time()
message()
message("Input annotation 1:             ", annot1_file)
message("Input annotation 2:             ", annot2_file)
message("Output annotation:              ", annot_out_file)
message("Output GO map file:             ", go_file)
message("Output KEGG map file:           ", kegg_file)
message("======================================================================")
message()


# CREATE A MERGED ANNOTATION FILE ----------------------------------------------
# Read the 1st annotation file
annot1 <- read_tsv(annot1_file,
                   col_types = cols(`EggNOG KEGG Terms` = col_character())) |>
  janitor::clean_names() |>
  # Remove genes with no sim.search hit -- those are also present in annot2
  filter(!is.na(subject_sequence))

# Read the 2nd annotation file (NR database only, with genes that had no sim.search hit in annot1)
annot2 <- read_tsv(annot2_file,
                   col_types = cols(`EggNOG KEGG Terms` = col_character())) |>
  janitor::clean_names()

# Check if there are no overlapping genes
stopifnot(!any(duplicated(c(annot1$query_sequence, annot2$query_sequence))))

# Row-bind the two annotations
annot <- bind_rows(annot1, annot2) |>
  select(!starts_with("x")) |>
  mutate(query_sequence = sub("t\\d+$", "", query_sequence),
         origin_database = case_when(
           grepl("refseq", origin_database, ignore.case = TRUE) ~ "refseq",
           grepl("uniprot", origin_database, ignore.case = TRUE) ~ "uniprot",
           grepl("_nr", origin_database, ignore.case = TRUE) ~ "nr"),
         ) |>
  rename(gene_id = query_sequence,
         best_hit = subject_sequence) |>
  arrange(gene_id) |>
  select(-contaminant)

# clean_names introduced some extra underscores we need to get rid of:
colnames(annot) <- sub("egg_nog", "eggnog", colnames(annot))
colnames(annot) <- sub("ip_scan", "ipscan", colnames(annot))
colnames(annot) <- sub("uni_prot", "uniprot", colnames(annot))


# PROCESS GO -------------------------------------------------------------------
progs <- c("eggnog", "uniprot", "ipscan")
cats <- c("biological", "cellular", "molecular")
go_cols <- paste0(rep(progs, each = length(cats)), "_go_", cats)

# Function to create GO-map for one column with GO terms
map_go <- function(col, annot) {
  # Get program and category
  program <- sub("_go_.*", "", col)
  ontology <- sub(".*_go_", "", col)
  
  # Get max. nr of terms
  n_terms <- max(str_count(annot[[col]], ",") + 1, na.rm = TRUE)

  annot |>
    select(gene_id, any_of(col), eggnog_tax_scope) |>
    separate_wider_delim(
      cols = any_of(col),
      delim = ",",
      too_few = "align_start",
      names = paste0("cat", seq_len(n_terms))
    ) |>
    pivot_longer(cols = starts_with("cat"), names_to = NULL, values_to = "category") |>
    filter(!is.na(category), category != "") |> 
    mutate(level = as.integer(sub(".*L=(-?\\d+)\\)$", "\\1", category)),
           category = sub("\\(L=.*", "", category),
           ontology = ontology,
           program = program,
           eggnog_tax_scope = ifelse(program == "eggnog", eggnog_tax_scope, NA),
           eggnog_tax_scope = as.character(eggnog_tax_scope)) |>
    relocate(category, .after = gene_id)
}

# Get all GO results in a single df with 1 gene-to-GOterm mapping per row
go_df <- map_dfr(.x = go_cols, .f = map_go, annot = annot) |>
  arrange(gene_id)


# PROCESS KEGG -------------------------------------------------------------------
n_terms <- max(str_count(annot$eggnog_kegg_terms, ",") + 1, na.rm = TRUE)

kegg_df <- annot |>
  select(gene_id, eggnog_kegg_terms, eggnog_tax_scope) |>
  filter(!is.na(eggnog_kegg_terms), eggnog_kegg_terms != "") |>
  separate_wider_delim(
    cols = eggnog_kegg_terms,
    delim = ",",
    too_few = "align_start",
    names = paste0("k", seq_len(n_terms))
  ) |>
  pivot_longer(cols = starts_with("k"), names_to = NULL, values_to = "category") |>
  filter(!is.na(category)) |>
  mutate(category = paste0("map", sprintf("%05d", as.integer(category)))) |>
  select(gene_id, category, eggnog_tax_scope) |>
  arrange(gene_id)


# GET AND REPORT SOME STATS ----------------------------------------------------
message("\n# Overview of database origin of similarity search annotations:")
table(annot$origin_database)

message("\n# Overview of taxonomic scope of EggNOG annotations:")
table(annot$eggnog_tax_scope)

n_tot <- nrow(annot)
n_egg <- sum(is.na(annot$best_hit))
n_inform <- sum(annot$informative == "Yes", na.rm = TRUE)
go_eggnog <- !is.na(annot$eggnog_go_biological) |
  !is.na(annot$eggnog_go_cellular) |
  !is.na(annot$eggnog_go_molecular)
go_uniprot <- !is.na(annot$uniprot_go_biological) |
  !is.na(annot$uniprot_go_cellular) |
  !is.na(annot$uniprot_go_molecular)
go_ipscan <- !is.na(annot$ipscan_go_biological) |
  !is.na(annot$ipscan_go_cellular) |
  !is.na(annot$ipscan_go_molecular)
n_go <- sum(go_eggnog | go_uniprot | go_ipscan)

kegg_uniprot <- !is.na(annot$uniprot_kegg_terms)
kegg_eggnog <- !is.na(annot$eggnog_kegg_terms)
n_kegg <- sum(kegg_eggnog | kegg_uniprot)

message()
message("======================================================================")
message("Total nr of genes:                                          ", n_tot)
message("Nr of genes with only EggNOG annotation:                    ", n_egg)
message("Nr of genes with informative similarity search annotation:  ", n_inform)
message()
message("Nr of genes with eggNOG GO terms:                           ", sum(go_eggnog))
message("Nr of genes with Uniprot GO terms:                          ", sum(go_uniprot))
message("Nr of genes with Interproscan GO terms:                     ", sum(go_ipscan))
message("Nr of genes with any GO term:                               ", n_go)
message("Number of gene-to-GOterm mapping:                           ", nrow(go_df))
message("Number of unique GO terms:                                  ", length(unique(go_df$category)))
message()
message("Nr of genes with eggNOG KEGG terms:                         ", sum(kegg_eggnog))
message("Nr of genes with Uniprot KEGG terms:                        ", sum(kegg_uniprot))
message("Nr of genes with any KEGG term:                             ", n_kegg)
message("Number of gene-to-KEGGterm mapping:                         ", nrow(kegg_df))
message("Number of unique KEGG terms:                                ", length(unique(kegg_df$category)))


# REMOVE GO AND KEGG COLUMNS FROM THE ANNOTATION -------------------------------
annot <- annot |> select(!contains("_go_")) |> select(!contains("_kegg_"))


# WRAP UP ----------------------------------------------------------------------
# Write output files
write_tsv(kegg_df, kegg_file)
write_tsv(go_df, go_file)
write_tsv(annot, annot_out_file)

# List output
message("\n# Listing the output file(s):")
system(paste("ls -lh", annot_out_file))
system(paste("ls -lh", go_file))
system(paste("ls -lh", kegg_file))

message("\n# Done with script entapnf_combine.R")
Sys.time()
message()
sessionInfo()
