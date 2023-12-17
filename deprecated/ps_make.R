#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=ps_make
#SBATCH --output=slurm-ps_make-%j.out

#? This script will build a phyloseq object, taking as input:
#? (1) A metadata table
#? (2) A table with counts for each sample and ASV ('seqtab' from 'dada.R')
#? (4) Taxonomic assignments for the ASVs in the dataset (from 'tax_assign_*.R')
#? (3) A phylogenetic tree of the ASVs in the dataset (from 'tree_build.R')
#? NOTE: In the metadata file, sample IDs should be in the 1st column

#? Load the Conda environment as follows to run this script directly using sbatch:
#? module load miniconda3/4.12.0-py39 && source activate /fs/ess/PAS0471/jelmer/conda/r-metabar

# SET-UP -----------------------------------------------------------------------
# Packages
packages <- c("BiocManager", "dada2", "phyloseq", "Biostrings", "QsRutils")

# Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--meta",
                    type = "character", required = TRUE,
                    help = "Input file with metadata (REQUIRED)")
parser$add_argument("--asv_table",
                    type = "character", required = TRUE,
                    help = "RDS file with sequence table (REQUIRED)")
parser$add_argument("--tree",
                    type = "character", required = TRUE,
                    help = "RDS file with tree (REQUIRED)")
parser$add_argument("--tax_table",
                    type = "character", required = TRUE,
                    help = "RDS file with tax. assignments (REQUIRED)")
parser$add_argument("-o", "--outfile",
                    type = "character", required = TRUE,
                    help = "Output phyloseq RDS file (REQUIRED)")
parser$add_argument("--sample_col_name",
                    type = "character", required = FALSE, default = "sampleID",
                    help = "Name that will be given to first column of the metadata, which should contain sample IDs")
args <- parser$parse_args()

seqtab_rds <- args$asv_table
taxa_rds <- args$tax_table
tree_rds <- args$tree
meta_file <- args$meta
ps_rds <- args$outfile
sample_col_name <- args$sample_col_name

# Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
if (!require(QsRutils)) remotes::install_github("QsRutils") 
pacman::p_load(char = packages)

# Create output dir if needed
outdir <- dirname(ps_rds)
dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)

# Report
message()
message("# ====================================================================")
message("#               STARTING SCRIPT PS_MAKE.R")
message("# ====================================================================")
Sys.time()
message("# ASV count table RDS file (input):   ", seqtab_rds)
message("# Taxa RDS file (input):              ", taxa_rds)
message("# Tree RDS file (input):              ", tree_rds)
message("# Sample data file (input):           ", meta_file)
message("# Phyloseq RDS file (output):         ", ps_rds)
message("# ====================================================================")
message()


# LOAD INPUT DATA --------------------------------------------------------------
# Sequence table from dada2
seqtab <- readRDS(seqtab_rds)
sampleIDs_seqtab <- sub("_S\\d+$", "", rownames(seqtab)) # Remove "_S16" suffix
rownames(seqtab) <- sampleIDs_seqtab

# Taxonomic assignments
taxa <- readRDS(taxa_rds)

# Tree
tree <- readRDS(tree_rds)

# Sample metadata
meta <- read.table(meta_file, sep = "\t", header = TRUE)
colnames(meta)[1] <- sample_col_name
rownames(meta) <- meta[[sample_col_name]]


# CHECK SAMPLE IDs -------------------------------------------------------------
# Check for matching sample names in FASTQ files and metadata file
message("\n# First 6 IDs from metadata:")
head(meta[[sample_col_name]])
message("# First 6 IDs from seqtab (i.e., from FASTQ file names):")
head(sampleIDs_seqtab)

message("\n# Are the sample IDs from the metadata and the seqtab the same?")
identical(sort(meta[[sample_col_name]]), sampleIDs_seqtab)

message("# Are any samples missing from the seqtab?")
setdiff(sort(meta[[sample_col_name]]), sampleIDs_seqtab)

message("# Are any samples missing from the metadata?")
setdiff(sampleIDs_seqtab, sort(meta[[sample_col_name]]))


# CREATE PHYLOSEQ OBJECT -------------------------------------------------------
ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               phy_tree(tree$tree),
               sample_data(meta),
               tax_table(taxa))

# Root the tree by the longest terminal branch (https://github.com/jfq3/QsRutils/blob/master/R/root_phyloseq_tree.R)
ps <- QsRutils::root_phyloseq_tree(ps)


# RENAME ASVs AND STORE SEQS SEPARATELY ----------------------------------------
# Extract ASV sequences (which are the taxa_names)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

# Merge sequence object into the phyloseq object:
ps <- merge_phyloseq(ps, dna)

# Rename ASVs from full sequences to ASV1...ASVx
taxa_names(ps) <- paste("ASV", 1:ntaxa(ps), sep = "_")


# SAVE AND WRAP UP -------------------------------------------------------------
# Save phyloseq object in an RDS file
saveRDS(ps, ps_rds)

# Report
message("\n# Listing output file:")
system(paste("ls -lh", ps_rds))
message("\n# Done with script ps_make.R")
Sys.time()
message()
sessionInfo()
message()
