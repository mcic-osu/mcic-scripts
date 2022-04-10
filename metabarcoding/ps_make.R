#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm-ps_make-%j.out

#? NOTE: In the metadata file, sample IDs should be in the 1st column

# SET-UP -----------------------------------------------------------------------
## Report
message()
message("## Starting script ps_make.R")
Sys.time()
message()

## Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-m", "--meta",
                    type = "character", required = TRUE,
                    help = "Input file with metadata (REQUIRED)")
parser$add_argument("-s", "--seqtab",
                    type = "character", required = TRUE,
                    help = "RDS file with sequence table (REQUIRED)")
parser$add_argument("-t", "--tree",
                    type = "character", required = TRUE,
                    help = "RDS file with tree (REQUIRED)")
parser$add_argument("-x", "--taxa",
                    type = "character", required = TRUE,
                    help = "RDS file with tax. assignments (REQUIRED)")
parser$add_argument("-o", "--outfile",
                    type = "character", required = TRUE,
                    help = "Output phyloseq RDS file (REQUIRED)")
args <- parser$parse_args()

seqtab_rds <- args$seqtab
taxa_rds <- args$taxa
tree_rds <- args$tree
meta_file <- args$meta
ps_rds <- args$outfile

## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "phyloseq", "DECIPHER")
pacman::p_load(char = packages)

## Create output dir if needed
outdir <- dirname(ps_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Report
message("## Sequence table RDS file (input):    ", seqtab_rds)
message("## Taxa RDS file (input):              ", taxa_rds)
message("## Tree RDS file (input):              ", tree_rds)
message("## Sample data file (input):           ", meta_file)
message("## Phyloseq RDS file (output):         ", ps_rds)
message()


# LOAD INPUT DATA --------------------------------------------------------------
## Sequence table from dada2
seqtab <- readRDS(seqtab_rds)
sampleIDs_seqtab <- sub("_S\\d+$", "", rownames(seqtab)) # Remove "_S16" suffix
rownames(seqtab) <- sampleIDs_seqtab

## Taxonomic assignments
taxa <- readRDS(taxa_rds)

## Tree
tree <- readRDS(tree_rds)

## Sample metadata
meta <- read.table(meta_file, sep = "\t", header = TRUE)
colnames(meta)[1] <- "sampleID"
rownames(meta) <- meta$sampleID


# CHECK SAMPLE IDs -------------------------------------------------------------
## Check for matching sample names in FASTQ files and metadata file
message("\n## First 6 IDs from metadata:")
head(meta$sampleID)
message("## First 6 IDs from seqtab (i.e., from FASTQ file names):")
head(sampleIDs_seqtab)

message("\n## Are the sample IDs from the metadata and the seqtab the same?")
identical(sort(meta$sampleID), sampleIDs_seqtab)

message("## Are any samples missing from the seqtab?")
setdiff(sort(meta$sampleID), sampleIDs_seqtab)

message("## Are any samples missing from the metadata?")
setdiff(sampleIDs_seqtab, sort(meta$sampleID))


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
message("\n## Listing output file:")
system(paste("ls -lh", ps_rds))
message("\n## Done with script ps_make.R")
Sys.time()
message()
