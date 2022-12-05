#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-tximport-%j.out

# DESCRIPTION ------------------------------------------------------------------
# This script will import Kallisto files into R and create a DEseq2 object
#TODO - Also accept Salmon counts
#TODO - Make column names for tx2gene/metadata explicit/flexible

# SET-UP -----------------------------------------------------------------------
message("\n## Starting script tximport.R")
Sys.time()
message()

## Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
if (!require(pacman)) install.packages("pacman", repos = rep, lib = lib)
if (!require(BiocManager)) install.packages("BiocManager", repos = rep, lib = lib)
if (!require(tximport)) BiocManager::install("tximport")
if (!require(rhdf5)) BiocManager::install("rhdf5")
if (!require(DESeq2)) BiocManager::install("DESeq2")

packages <- c("argparse", "tximport", "rhdf5", "DESeq2", "tidyverse")
pacman::p_load(char = packages, install = TRUE)

## Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--indir",
                    type = "character", required = TRUE, default = NULL,
                    help = "Input dir with Kallisto output files (REQUIRED)")
parser$add_argument("-o", "--outdir",
                    type = "character", required = TRUE, default = NULL,
                    help = "Output directory")
parser$add_argument("-m", "--meta",
                    type = "character", required = TRUE, default = NULL,
                    help = "Metadata file (TSV)")
parser$add_argument("-t", "--tx2gene",
                    type = "integer", required = TRUE, default = NULL,
                    help = "Transcript-to-gene map (TSV)")
args <- parser$parse_args()

## Output file
dds_out <- here(outdir, "deseq_object.rds")

## Settings
dirname_column <- "id_long"    # Column in metadata with Kallisto dirnames
sampleid_column <- "id_short"  # Column in metadata with desired sample names


# PROCESS METADATA AND KALLISTO COUNTS -----------------------------------------
## Read metadata
meta <- read.delim(meta_file) %>% arrange(.data[[sampleid_column]])
rownames(meta) <- meta[[sampleid_column]]

## Read transcript-to-gene map
transmap <- read_tsv(transmap_file,
                     col_names = c("GENEID", "TXNAME"),
                     col_types = "cc") %>%
  select(TXNAME, GENEID)

## Import Kallisto transcript counts --
## create gene-level count estimates normalized by library size and transcript length
## See https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
kallisto_files <- here(kallisto_dir, meta[[dirname_column]], "abundance.h5")
names(kallisto_files) <- meta[[sampleid_column]]

txi <- tximport(kallisto_files,
                type = "kallisto",
                tx2gene = transmap,
                countsFromAbundance = "lengthScaledTPM")


# CREATE DESEQ OBJECT ----------------------------------------------------------
## Check that sample names are the same and samples are in same order
stopifnot(all(rownames(meta) == colnames(txi$counts)))
message("\n## Sample names:")
print(rownames(meta))
message("\n## Dimensions of count matrix:")
dim(txi$counts)

## Create DESeq object
dds <- DESeqDataSetFromTximport(txi, meta, ~1)

## Save DESeq object
saveRDS(dds, dds_out)


# WRAP UP ----------------------------------------------------------------------
## List output
message("\n## Listing output file:")
system(paste("ls -lh", dds_out))

message("\n## Done with script tximport.R")
Sys.time()
message()
