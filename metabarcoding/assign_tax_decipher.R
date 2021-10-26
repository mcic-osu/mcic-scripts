# SET-UP -----------------------------------------------------------------------
## Load packages
if(!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("dada2", "DECIPHER")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
seqtab_rds <- args[1]
taxa_rds <- args[2]
n_cores <- as.integer(args[3])

# seqtab_rds <- "results/ASV/main/seqtab.rds"
# taxa_rds <- "results/ASV/main/taxa_decipher.rds"
# n_cores <- 8

## Create output dir if needed
outdir <- dirname(taxa_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Other variables/constants
tax_URL <- "http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData"
tax_file <- file.path(outdir, basename(tax_URL))
species_URL <- "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"
species_file <- file.path(outdir, basename(species_URL))

qc_file <- file.path(outdir, "tax_prop_assigned_decipher.txt")

## Report
cat("## Starting script assign_tax_decipher.R\n")
Sys.time()
cat("## Sequence table RDS file (input):", seqtab_rds, "\n")
cat("## Taxa RDS file (output):", taxa_rds, "\n")
cat("## Proportion-assigned QC file (output):", qc_file, "\n")
cat("## Number of cores:", n_cores, "\n\n")
cat("## Taxonomic assignment file (downloaded input):", tax_file, "\n")
cat("## Species assignment file (downloaded input):", species_file, "\n\n")

## Function to get the proportion of ASVs assigned to taxa
qc_tax <- function(taxa) {
    prop <- apply(taxa, 2,
                  function(x) round(length(which(!is.na(x))) / nrow(taxa), 4))
    prop <- data.frame(prop)
    prop$tax_level <- rownames(prop) 
    rownames(prop) <- NULL
    return(prop)
}

# PREPARE INPUT DATA -----------------------------------------------------------
## Create a DNAStringSet from the ASVs
seqtab <- readRDS(seqtab_rds)
dna <- DNAStringSet(getSequences(seqtab))


# DECIPHER TAXONOMIC ASSIGNMENT ------------------------------------------------
cat("## Assigning taxonomic labels to ASVs...\n")

## Get and load DECIPHER training set
### https://www.bioconductor.org/packages/devel/bioc/vignettes/DECIPHER/inst/doc/ClassifySequences.pdf
### http://www2.decipher.codes/Downloads.html
if (!file.exists(tax_file)) download.file(url = tax_URL, destfile = tax_file)
load(tax_file)   # Will create object "trainingSet"
if (!file.exists(species_file)) download.file(url = species_URL, destfile = species_file)

## Assign taxonomy
cat("\n## Now assigning taxonomy...\n")
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

## addSpecies doesn't work after idtaxa!
#cat("\n## Now adding species-level assignments...\n")
#taxa <- addSpecies(taxa, species_file)

## QC: Check proportions of ASVs assigned to taxa
prop_assigned <- qc_tax(taxa)
cat("\n## Proportion of ASVs assigned to different taxonomic levels - DECIPHER:\n")
print(prop_assigned)
write.table(prop_assigned, qc_file, sep = "\t", quote = FALSE, row.names = FALSE)

## Save RDS file
saveRDS(taxa, taxa_rds)


# WRAP UP ----------------------------------------------------------------------
## Report
cat("\n## Done with script assign_tax_decipher.R\n")
cat("## Taxa RDS output file:", taxa_rds, "\n")
cat("## Proportion-assigned QC file:", qc_file, "\n")
Sys.time()
