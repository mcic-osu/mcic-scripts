# SET-UP -----------------------------------------------------------------------
## Load packages
if(!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("dada2", "phyloseq", "DECIPHER")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
seqtab_rds <- args[1]
taxa_rds <- args[2]
tree_rds <- args[3]
sampledata_file <- args[4]
ps_rds <- args[5]

## Example parameters for interactive testing:
# seqtab_rds <- "results/ASV/main/seqtab.rds"
# taxa_rds <- "results/taxonomy/taxa_dada.rds"
# tree_rds <- "results/ASV/main/tree.rds"
# sampledata_file <- "metadata/sample_data.txt"
# ps_rds <- "results/ASV/main/ps.rds"

## Create output dir if needed
outdir <- dirname(ps_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Report
cat("## Starting script make_ps.R\n")
Sys.time()
cat("## Sequence table RDS file (input):", seqtab_rds, "\n")
cat("## Taxa RDS file (input):", taxa_rds, "\n")
cat("## Tree RDS file (input):", tree_rds, "\n")
cat("## Sample data file (input):", sampledata_file, "\n")
cat("## Phyloseq RDS file (output):", ps_rds, "\n\n")

# LOAD INPUT DATA --------------------------------------------------------------
## Sequence table from dada2
seqtab <- readRDS(seqtab_rds)
sampleIDs_seqtab <- sub("_S\\d+$", "", rownames(seqtab)) # Remove "_S16" etc suffix
rownames(seqtab) <- sampleIDs_seqtab

## Taxonomic assignments
taxa <- readRDS(taxa_rds) 

## Tree
tree <- readRDS(tree_rds)

## Sample metadata
meta <- read.table(sampledata_file, sep = "\t", header = TRUE)
colnames(meta)[1] <- "sample_ID"
rownames(meta) <- meta$sample_ID

# CHECK SAMPLE IDs -------------------------------------------------------------
## Check for matching sample names in FASTQ files and metadata file
cat("## IDs from metadata:\n")
head(meta$sample_ID)
cat("## IDs from seqtab (i.e., from FASTQ file names):\n")
head(sampleIDs_seqtab)

cat("\n## Are the sample IDs from the metadata and the seqtab the same?\n")
identical(sort(meta$sample_ID), sampleIDs_seqtab)

cat("## Are any samples missing from the seqtab?\n")
setdiff(sort(meta$sample_ID), sampleIDs_seqtab)

cat("## Are any samples missing from the metadata?\n")
setdiff(sampleIDs_seqtab, sort(meta$sample_ID))


# CREATE PHYLOSEQ OBJECT -------------------------------------------------------
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               phy_tree(tree$tree),
               sample_data(meta),
               tax_table(taxa))


# RENAME ASVs AND STORE SEQS SEPARATELY ----------------------------------------
## Extract ASV sequences (which are the taxa_names)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

## Merge sequence object into the phyloseq object:
ps <- merge_phyloseq(ps, dna)

## Rename ASVs from full sequences to ASV1...ASVx
taxa_names(ps) <- paste("ASV", 1:ntaxa(ps), sep = "_")


# SAVE AND WRAP UP -------------------------------------------------------------
## Save phyloseq object in an RDS file
saveRDS(ps, ps_rds)

## Report
q
