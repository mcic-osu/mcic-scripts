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
#? (2) A table with counts for each sample and ASV
#? (4) Taxonomic assignments for the ASVs in the dataset
#? (3) A phylogenetic tree of the ASVs in the dataset
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
args <- parser$parse_args()

asv_table_file <- args$asv_table
tax_table_file <- args$tax_table
tree_file <- args$tree
meta_file <- args$meta
ps_rds <- args$outfile

#asv_table_file <- "16S/results/ampliseq_16S/qiime2/abundance_tables/feature-table.tsv"
#tax_table_file <- "16S/results/ampliseq_16S/dada2/ASV_tax_species.tsv"
#tree_file <- "16S/results/ampliseq_16S/qiime2/phylogenetic_tree/tree.nwk"
#meta_file <- "16S/data/meta/sample_data.tsv"
#ps_rds <-"results/phyloseq/ps.rds"

# Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
if (!require(QsRutils)) remotes::install_github("jfq3/QsRutils")
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
message("# ASV count table RDS file (input):   ", asv_table_file)
message("# Taxa RDS file (input):              ", tax_table_file)
message("# Tree RDS file (input):              ", tree_file)
message("# Sample data file (input):           ", meta_file)
message("# Phyloseq RDS file (output):         ", ps_rds)
message("# ====================================================================")
message()


# LOAD INPUT DATA --------------------------------------------------------------
# Sequence table from dada2
asv_mat <- read.delim(asv_table_file, skip=1, comment.char = "", sep = "\t")
rownames(asv_mat) <- asv_mat[[1]]
asv_mat <- asv_mat[, -1]

# Taxonomic assignments
tax_mat <- read.delim(tax_table_file, sep = "\t")
rownames(tax_mat) <- tax_mat[[1]]
tax_mat <- tax_mat[, -1]
tax_mat <- as.matrix(tax_mat)
sequences <- tax_mat[, "sequence"] # Put sequences in separate table
tax_mat <- tax_mat[, -which(colnames(tax_mat) == "sequence")] # Get rid of sequence
tax_mat <- gsub("^$", NA, tax_mat) # Empty strings should be NA

# Tree
tree <- read_tree(tree_file)

# Sample metadata
meta <- read.table(meta_file, sep = "\t", header = TRUE)
rownames(meta) <- meta[[1]]
meta <- meta[, -1]


# CREATE PHYLOSEQ OBJECT -------------------------------------------------------
ps <- phyloseq(otu_table(asv_mat, taxa_are_rows = TRUE),
               tax_table(tax_mat),
               phy_tree(tree),
               sample_data(meta),
               Biostrings::DNAStringSet(sequences))

# Root the tree by the longest terminal branch (https://github.com/jfq3/QsRutils/blob/master/R/root_phyloseq_tree.R)
ps <- QsRutils::root_phyloseq_tree(ps)


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
