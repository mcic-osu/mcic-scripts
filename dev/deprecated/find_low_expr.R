#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=filter_expr
#SBATCH --output=slurm-filter_expr-%j.out

# Parse command-line args
args <- commandArgs(trailingOnly = TRUE)
mat_file <- args[1] #mat_file <- "/fs/ess/PAS0471/jelmer/assist/2022-04_soumya/TRAM_SUBSET/results/filter_expr/norm_alltrans.gene.TPM.not_cross_norm"
outfile <- args[2]
max_tpm <- as.numeric(args[3])
mean_tpm <- as.numeric(args[4])

# Report
message("# Starting script find_low_expr.R")
message("Count matrix file:                     ", mat_file)
message("Output file:                           ", outfile)
message("Min-max TPM:                           ", max_tpm)
message("Min-mean TPM:                          ", mean_tpm, "\n")

# Read the counts file
counts <- read.delim(mat_file, row.names = 1)

# Check which genes should be filtered
max_too_low <- names(which(apply(counts, 1, max) < 1))
mean_too_low <- names(which(rowMeans(counts) < 0.1))
either_too_low <- unique(c(max_too_low, mean_too_low))

# Write the output file
writeLines(either_too_low, outfile)

# Report
message("# Statistics:")
message("Number of input genes:                 ", nrow(counts))
message("Number of genes with max TPM<1:        ", length(max_too_low))
message("Number of genes with mean TPM<0.1:     ", length(mean_too_low))
message("Number of genes to be removed:         ", length(either_too_low))
message("Proportion of genes to be removed:     ", round(length(either_too_low) / nrow(counts), 3))
message("Number of genes to be kept:            ", nrow(counts) - length(either_too_low))
message()
message("Done with script")
Sys.time()
message()
