#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-assigntax4qiime-%j.out

# SETUP ------------------------------------------------------------------------
## Load packages
if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org/")
packages <- c("tidyverse",
              "qiime2R",      # https://github.com/jbisanz/qiime2R
              "dada2",
              "argparse")
pacman::p_load(char = packages, install = TRUE, repos = "https://cloud.r-project.org/")

## Parse options
parser <- ArgumentParser() # create parser object
parser$add_argument("-i", "--seq_artifact", type = "character", default = NULL, 
                    help = "Input tree file (REQUIRED)")
parser$add_argument("-o", "--outdir", type = "character", default = NULL, 
                    help = "Output dir (REQUIRED)")
parser$add_argument("-l", "--locus", type = "character", default = NULL,
                    help = "Locus type: '16S' or 'ITS' (REQUIRED)")
parser$add_argument("-r", "--ref_file", type = "character", default = NULL, 
                    help = "FASTA file with reference sequences from UNITE/SILVA/etc")
parser$add_argument("-s", "--species_file", type = "character", default = NULL, 
                    help = "Silva species file (16S only)")
args <- parser$parse_args()

seq_artifact <- args$seq_artifact
outdir <- args$outdir
locus <- args$locus
ref_file <- args$ref_file
species_file <- args$ref_file

## Constants
TAX_LEVELS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

## Define output files
tax_rds <- file.path(outdir, "tax.rds")
tax_raw_rds <- file.path(outdir, "tax_rawoutput.rds")
tax_txt <- file.path(outdir, "tax.txt")
qc_file <- file.path(outdir, "prop_assigned.txt")
plot_file <- file.path(outdir, "prop_assigned.png")

## Create output dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## Report
message("## Starting script assigntax4qiime.R...")
Sys.time()
message()
message("## Qiime2 sequence artifact:            ", seq_artifact)
message("## Output dir:                          ", outdir)
message("## Locus type (16S or ITS):             ", locus)
if (!is.null(ref_file)) message("## Reference FASTA:                     ", ref_file)
if (!is.null(species_file)) message("## Species reference FASTA:             ", species_file)
message("--------------------\n")

# PREP REF FILES ---------------------------------------------------------------
if (locus == "16S") {
  if (is.null(ref_file)) {
    ref_url <- "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
    ref_file <- file.path(outdir, basename(ref_url))
    message("## Taxonomic reference file:            ", ref_file)
    if (!file.exists(ref_file)) {
      message("## Downloading Silva reference sequences...")
      download.file(url = ref_url, destfile = ref_file)
    }
  }
  if (is.null(species_file)) {
    species_url <- "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"
    species_file <- file.path(outdir, basename(species_url))
    message("## Species reference file:              ", species_file)
    if (!file.exists(species_file)) {
      message("## Downloading Silva reference sequences for species assignment...")
      download.file(url = species_url, destfile = species_file)
    }
  }
}


# MAIN -------------------------------------------------------------------------
## Read input files
seqs <- read_qza(seq_artifact)$data
seqs_df <- data.frame(qiimeID = names(seqs), seq = as.character(seqs))

## Assign taxonomy
message("\n## Now assigning taxonomy...")
taxa <- assignTaxonomy(seqs, ref_file, tryRC = TRUE, multithread = 8)
if (locus == "16S") {
  message("\n## Adding species-level taxonomy...")
  taxa <- addSpecies(taxa, species_file)
}
saveRDS(taxa, tax_raw_rds)


# POST-PROCESS -----------------------------------------------------------------
message("\n## Now processing assignment output...")
if (locus == "ITS") {
  ## Process output for usage in R - get rid of "p__", "k__", etc prefixes
  tax_df <- taxa %>%
    as.data.frame() %>%
    mutate(across(everything(), ~ sub("[kpcofgs]__", "", .x)))
  saveRDS(tax_df, tax_rds)
}

## Process for Qiime2 import - https://forum.qiime2.org/t/importing-taxonomy-tables-from-dada2/3609/10
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- as.data.frame(taxa) %>%
  unite(col = "taxonomy", all_of(tax_cols), sep = ";") %>%
  merge(seqs_df, by.x = "row.names", by.y = "seq") %>%
  select(qiimeID, taxonomy)

## Write to text file that can be imported by Qiime
write.table(tax, tax_txt,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


# QC ---------------------------------------------------------------------------
## Function to get the proportion of ASVs assigned to taxa
qc_tax <- function(taxa,
                   tax_levels = c("Kingdom", "Phylum", "Class", "Order", "Family",
                                  "Genus", "Species")) {
  n <- apply(taxa, 2, function(x) length(which(!is.na(x))))
  prop <- round(n / nrow(taxa), 4)
  data.frame(n, prop) %>%
    rownames_to_column("tax_level") %>%
    mutate(tax_level = factor(tax_level, levels = tax_levels))
}

## Create df with proportion assigned
prop_assigned <- qc_tax(taxa, tax_levels = TAX_LEVELS)
write_tsv(prop_assigned, qc_file)
message("\n## Prop. ASVs assigned to taxonomy:")
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
message("\n## Listing output files:")
system(paste("ls -lh", outdir))
message("\n## Done with script.")
Sys.time()
