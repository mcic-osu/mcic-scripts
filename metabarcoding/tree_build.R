# SET-UP -----------------------------------------------------------------------
## Load packages
if(!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("dada2", "DECIPHER", "phangorn")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
seqtab_rds <- args[1]             # Sequence table RDS file (input)
tree_rds <- args[2]               # Tree RDS file (output)
n_cores <- as.integer(args[3])    # Number of computer cores to use

# seqtab_rds <- "results/dada2/main/rds/seqtab_nonchim_lenfilter.rds"
# tree_rds <- "results/tree/tree.rds"

## Report
cat("## Starting script tree.R\n")
Sys.time()
cat("## Sequence table RDS file:", seqtab_rds, "\n")
cat("## Tree RDS file:", tree_rds, "\n")
cat("## Number of cores:", n_cores, "\n\n")

## Create output dir if needed
outdir <- dirname(tree_rds)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Load the input data
seqtab <- readRDS(seqtab_rds)

seqs <- getSequences(seqtab)
names(seqs) <- seqs


# BUILD THE TREE ---------------------------------------------------------------
# See https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.html

## 1 - Align
cat("## Aligning sequences...\n")
alignment <- AlignSeqs(DNAStringSet(seqs),
                       anchor = NA,
                       iterations = 5,
                       refinements = 5,
                       processors = n_cores)

## 2 - Compute distances
cat("## Computing pairwise distances from ASVs...\n")
alignment_mat <- phyDat(as(alignment, "matrix"), type = "DNA")
dist_mat <- dist.ml(alignment_mat)

## 3 - Build neighbor-joining tree and compute its likelihood
cat("## Building a tree...\n")
treeNJ <- NJ(dist_mat)                     # Build NJ tree
fit <- pml(treeNJ, data = alignment_mat)   # Compute likelihood
fitGTR <- update(fit, k = 4, inv = 0.2)    # Update to GTR model

## 4 - Compute likelihood
cat("## Optimizing the tree...\n")
fitGTR <- optim.pml(fitGTR,
                    model = "GTR",
                    optInv = TRUE,
                    optGamma = TRUE,
                    rearrangement = "stochastic",
                    control = pml.control(trace = 0))


# WRAP UP ----------------------------------------------------------------------
## Save RDS file
saveRDS(fitGTR, tree_rds)

## Report
cat("\n## Listing output files:\n")
system(paste("ls -lh", tree_rds))

cat("## Done with script tree.R\n")
Sys.time()
