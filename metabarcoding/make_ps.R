# SET-UP --------------------------------------

## Set number of cores
n_cores <- 8

## Load packages
print("Loading packages...")
if(!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse", "gridExtra", "dada2",
              "phyloseq", "DECIPHER", "phangorn")
pacman::p_load(char = packages)

## Define dirs
refdata_dir <- "data/ref"                # For reference data like tax db's

if (!dir.exists(refdata_dir)) dir.create(refdata_dir, recursive = TRUE)

## Define FASTA file with training data
## (Check for an up-to-date version at <https://benjjneb.github.io/dada2/training.html>)
tax_file <- file.path(refdata_dir, "silva_nr99_v138.1_train_set.fa.gz")
tax_URL <- "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
if (!file.exists(tax_file)) download.file(url = tax_URL, destfile = tax_file)

## Define file with sample metadata
metadata_file <- "metadata/sample_data.txt"


# PREP METADATA ----------------------------------------------------

## Load sample metadata
print("Load and prepare sample metadata...")

metadata_df <- read.table(file = metadata_file, sep = "\t", header = TRUE)
colnames(metadata_df)[1] <- "SampleID"
rownames(metadata_df) <- metadata_df$SampleID

print("Head of metadata file:")
head(metadata_df)

## Check for matching sample names in FASTQ files and metadata file
print("IDs from metadata:")
metadata_df$SampleID
print("FASTQ file names:")
head(basename(fastqs_raw_F))

print("Are the sample IDs from the metadata and the fastq files the same?")
identical(sort(metadata_df$SampleID), sampleIDs)

print("Are any samples missing from the fastq files?")
setdiff(sort(metadata_df$SampleID), sampleIDs)

print("Are any samples missing from the metadata?")
setdiff(sampleIDs, sort(metadata_df$SampleID))


# ASSIGN TAXONOMY ------------------------------------------------

print("Assigning taxonomic labels to ASVs...")

taxa <- assignTaxonomy(seqtab, tax_key, multithread = n_cores)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

saveRDS(taxa, file = file.path(outdir, "taxa.rds"))
# taxa <- readRDS(file.path(outdir, "taxa.rds"))

# CREATE PHYLOSEQ OBJECT -------------------------------------------

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(metadata_df),
               tax_table(taxa))

saveRDS(ps, file = file.path(outdir, "ps_V4.rds"))


print("Done with script.")