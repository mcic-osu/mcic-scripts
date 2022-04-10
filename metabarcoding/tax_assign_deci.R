#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-tax-assign-dada-%j.out


# SET-UP -----------------------------------------------------------------------
## Report
message("\n## Starting script tax_assign_deci.R")
Sys.time()
message()

## Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--seqtab_in",
                    type = "character", required = TRUE,
                    help = "Input file (sequence table RDS) (REQUIRED)")
parser$add_argument("-o", "--taxa_out",
                    type = "character", required = TRUE,
                    help = "Output file (taxa RDS file) (REQUIRED)")
parser$add_argument("-r", "--ref_url",
                    type = "character",
                    default = "http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData",
                    help = "Taxonomic reference URL [default %(default)s]")
parser$add_argument("-c", "--cores",
                    type = "integer", default = 1,
                    help = "Number of cores (threads) to use [default %(default)s]")
args <- parser$parse_args()

seqtab_rds <- args$seqtab_in
taxa_rds <- args$taxa_out
ref_url <- args$ref_url
n_cores <- args$cores

## Load other packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
packages <- c("BiocManager", "dada2", "DECIPHER", "tidyverse")
pacman::p_load(char = packages)

## Constants
TAX_LEVELS <- c("domain", "phylum", "class", "order", "family", "genus", "species")

## Create output dir if needed
outdir <- dirname(taxa_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Define output files
tax_file <- file.path(outdir, basename(ref_url))
qc_file <- file.path(outdir, "tax_prop_assigned_decipher.txt")
plot_file <- file.path(outdir, "tax_prop_assigned_decipher.png")

## Report
message("## Sequence table RDS file (input):               ", seqtab_rds)
message("## Taxa RDS file (output):                        ", taxa_rds)
message("## Proportion-assigned QC file (output):          ", qc_file)
message("## Number of cores:                               ", n_cores)
message()
message("## Taxonomic assignment file (downloaded input):  ", tax_file)
message()


# FUNCTIONS --------------------------------------------------------------------
## Function to get the proportion of ASVs assigned to taxa
qc_tax <- function(taxa, TAX_LEVELS) {
    prop <- apply(taxa, 2,
                  function(x) round(length(which(!is.na(x))) / nrow(taxa), 4))
    
    prop <- data.frame(prop) %>%
        rownames_to_column("tax_level") %>%
        mutate(tax_level = factor(tax_level, levels = TAX_LEVELS))

    return(prop)
}


# PREPARE INPUT DATA -----------------------------------------------------------
## Create a DNAStringSet from the ASVs
seqtab <- readRDS(seqtab_rds)
dna <- DNAStringSet(getSequences(seqtab))


# DECIPHER TAXONOMIC ASSIGNMENT ------------------------------------------------

## Get and load DECIPHER training set
#? https://www.bioconductor.org/packages/devel/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
#? http://www2.decipher.codes/Downloads.html
if (!file.exists(tax_file)) download.file(url = ref_url, destfile = tax_file)
load(tax_file)   # Will create object "trainingSet"

## Assign taxonomy
message("\n## Assigning taxonomy...")
ids <- IdTaxa(dna, trainingSet,
              strand = "top",
              threshold = 60,
              processors = n_cores, verbose = TRUE)

## Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab)
taxa <- taxid

## Save RDS file
saveRDS(taxa, taxa_rds)


# QC TAX. ASSIGNMENTS ----------------------------------------------------------
## Create df with proportion assigned
prop_assigned <- qc_tax(taxa, TAX_LEVELS = TAX_LEVELS)
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
system(paste("ls -lh", taxa_rds))
system(paste("ls -lh", qc_file))
system(paste("ls -lh", plot_file))

message("\n## Done with script tax_assign_deci.R")
Sys.time()
message()
