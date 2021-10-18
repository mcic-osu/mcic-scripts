#!/usr/bin/env Rscript

# SET-UP --------------------------------------

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)

fastq_indir <- args[1]          # Dir with input FASTQ files
outdir <- args[2]               # Dir for output 
n_cores <- as.integer(args[3])  # Number of computer cores to use
n_samples <- args[4]            # Number of samples to run ("all" => all samples, otherwise specify an integer)
trunc_f <- as.integer(args[5])  # Truncate F reads after trunc_f bases
trunc_r <- as.integer(args[6])  # Truncate R reads after trunc_r bases

# fastq_indir <- "results/cutadapt"
# outdir <- "results/ASV_inference"
# n_cores <- 8
# n_samples <- "5"
# trunc_f <- 180
# trunc_r <- 180

## Constants
save_rds <- TRUE               # Whether to save RDS files after every step

## Report command-line arguments
cat("## Starting script ASV_inference.R...")
Sys.time()
cat("## Dir with input FASTQ files:", fastq_indir, "\n")
cat("## Output dir:", outdir, "\n")
cat("## Number of cores:", n_cores, "\n")
cat("## Number of samples to analyze:", n_samples, "\n")
cat("## Truncate F after n bases:", trunc_f, "\n")
cat("## Truncate R after n bases:", trunc_r, "\n\n")

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse", "gridExtra", "dada2",
              "phyloseq", "DECIPHER", "phangorn")
pacman::p_load(char = packages)

## Define and create output dirs
filter_dir <- file.path(outdir, "fastq_filtered")   # For filtered FASTQ files

if (!dir.exists(filter_dir)) dir.create(filter_dir, recursive = TRUE)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Get paths to input FASTQ files
fastqs_raw_F <- sort(list.files(fastq_indir, pattern = "_R1_001.fastq.gz",
                                full.names = TRUE))
fastqs_raw_R <- sort(list.files(fastq_indir, pattern = "_R2_001.fastq.gz",
                                full.names = TRUE))

n_samples_all <- length(fastqs_raw_F)

if (n_samples == "all") {
  n_samples <- n_samples_all
} else {
  n_samples <- as.integer(n_samples)
}

fastqs_raw_F <- fastqs_raw_F[1:n_samples]
fastqs_raw_R <- fastqs_raw_R[1:n_samples]

cat("## Number of samples to be analyzed:", n_samples, "\n")
cat("## First 6 FASTA files:", head(fastqs_raw_F), "\n")

## Extract sample IDs from FASTQ file names 
sampleIDs <- sub("_L001.*", "", basename(fastqs_raw_F))

## Define output files
fastqs_filt_F <- file.path(filter_dir, paste0(sampleIDs, "_F_filt.fastq"))
fastqs_filt_R <- file.path(filter_dir, paste0(sampleIDs, "_R_filt.fastq"))

errorplot_F_file <- file.path(outdir, "errors_F.png")
errorplot_R_file <- file.path(outdir, "errors_R.png")

fasta_out <- file.path(outdir, "ASVs.fa")
seqtab_file <- file.path(outdir, "seqtab.rds")
qc_file <- file.path(outdir, "nreads_summary.txt")


# FASTQ FILE FILTERING  --------------------------------------------
cat("\n----------------\n## Step 1: Filtering and trimming FASTQ files...\n")

filter_results <-
  filterAndTrim(fastqs_raw_F, fastqs_filt_F,
                fastqs_raw_R, fastqs_filt_R,
                truncLen = c(trunc_f, trunc_r), #c(180, 160), 
                trimLeft = 0,
                trimRight = 0,
                maxN = 0,
                maxEE = c(2,2),
                truncQ = 2,
                rm.phix = FALSE,
                multithread = n_cores, 
                compress = FALSE, verbose = TRUE) 

head(filter_results)


# FASTQ FILE DEREPLICATION ------------------------------------------
cat("\n----------------\n## Step 2: Dereplicating FASTQ files...\n")

fastqs_derep_F <- derepFastq(fastqs_filt_F, verbose = FALSE)
fastqs_derep_R <- derepFastq(fastqs_filt_R, verbose = FALSE)

names(fastqs_derep_F) <- sampleIDs
names(fastqs_derep_R) <- sampleIDs

## Save objects to RDS files
if (save_rds) saveRDS(fastqs_derep_F, file = file.path(outdir, "fastqs_derep_F.rds"))
if (save_rds) saveRDS(fastqs_derep_R, file = file.path(outdir, "fastqs_derep_R.rds"))


# ERROR LEARNING ----------------------------------------------------
cat("\n----------------\n## Step 3: Learning errors...\n")

err_F <- learnErrors(fastqs_derep_F, multithread = n_cores, verbose = TRUE)
err_R <- learnErrors(fastqs_derep_R, multithread = n_cores, verbose = TRUE)

## Save objects to RDS files
if (save_rds) saveRDS(err_F, file = file.path(outdir, "err_F.rds"))
if (save_rds) saveRDS(err_R, file = file.path(outdir, "err_R.rds"))

## Plot errors
p <- plotErrors(err_F, nominalQ = TRUE)
ggsave(errorplot_F_file, width = 8, height = 7)
p <- plotErrors(err_R, nominalQ = TRUE)
ggsave(errorplot_R_file, width = 8, height = 7)


# INFER ASVS ------------------------------------------------ 
cat("\n----------------\n## Step 4: Inferring ASVs...\n")

dada_Fs <- dada(fastqs_derep_F, err = err_F, pool = FALSE, multithread = n_cores)
dada_Rs <- dada(fastqs_derep_R, err = err_R, pool = FALSE, multithread = n_cores)

## TODO - USE/TRY POOL=TRUE

## Save objects to RDS files
if (save_rds) saveRDS(dada_Fs, file = file.path(outdir, "dada_Fs.rds"))
if (save_rds) saveRDS(dada_Rs, file = file.path(outdir, "dada_Rs.rds"))


# MERGE READ PAIRS -----------------------------------------
cat("\n----------------\n## Step 5: Merging read pairs...\n")

mergers <- mergePairs(dada_Fs, fastqs_derep_F,
                      dada_Rs, fastqs_derep_R,
                      verbose = TRUE)

## Save objects to RDS files
if (save_rds) saveRDS(mergers, file = file.path(outdir, "mergers.rds"))


# CREATE SEQUENCE TABLE --------------------------------------------------
cat("\n----------------\n## Step 6: Creating the sequence table...\n")

seqtab_all <- makeSequenceTable(mergers)

## The dimensions of the object are the nr of samples (rows) and the nr of ASVs (columns):
cat("## Nr of ASVs before chimera removal:", ncol(seqtab_all), "\n")


# REMOVE CHIMERAS ----------------------------------------------
cat("\n----------------\n## Step 7: Removing chimeras...\n")

seqtab <- removeBimeraDenovo(seqtab_all,
                             method = "consensus",
                             multithread = n_cores,
                             verbose = TRUE)

cat("## Nr of ASVs after chimera removal:", ncol(seqtab), "\n")

## Save objects to RDS files
if (save_rds) saveRDS(seqtab, file = seqtab_file)


# CHECK ASV SEQUENCE LENTGTHS ---------------------------------
cat("\n----------------\n## Table of sequence lengths:\n")
table(nchar(getSequences(seqtab)))

## If you need to remove sequences of a particular length (e.g. too long):
## seqtab2 <- seqtab[, nchar(colnames(seqtab_all)) %in% seq(250,256)]

## TODO - calculate abundance sums for ASVs of different lengths

# CREATE QC SUMMARY TABLE -----------------------------------
cat("\n----------------\n## Creating QC summary table...\n")

## Define a function to get the nr of unique sequences
getN <- function(x) sum(getUniques(x))

## Calculate nr of unique sequences across denoised and merged seqs
denoised_F <- sapply(dada_Fs, getN)
denoised_R <- sapply(dada_Rs, getN)
merged <- sapply(mergers, getN)

## Put together the final QC table
nreads_summary <- data.frame(filter_results,
                             denoised_F,
                             denoised_R,
                             merged,
                             nonchim = rowSums(seqtab),
                             row.names = sampleIDs)
colnames(nreads_summary)[1:2] <- c("input", "filtered")

## Have a look at the first few rows
cat("## First few rows of the QC table:\n")
head(nreads_summary)

## Write QC summary table to file
write.table(nreads_summary,
            file = qc_file,
            sep = "\t", quote = FALSE, row.names = TRUE)


# CREATE AND WRITE FASTA FILE ----------------------------------------------
cat("\n----------------\n## Writing FASTA file...\n")

## Prepare sequences and headers
asv_seqs <- colnames(seqtab)
asv_headers <- paste(">ASV", 1:ncol(seqtab), sep = "_")

## Interleave headers and sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))

## Write FASTA file
write(asv_fasta, file = fasta_out)


# LIST OUTPUT FILES --------------------------------------------------------
cat("\n----------------\n## Listing output files:\n")

cat("## First few filtered FASTQ files:\n", fastqs_filt_F[1:2], "\n")

cat("## Error profile plot - F:", errorplot_F_file, "\n")
cat("## Error profile plot - R:", errorplot_R_file, "\n")

cat("## FASTA:", fasta_out, "\n")
cat("## Sequence table:", seqtab_file, "\n")
cat("## QC table:", qc_file, "\n")

cat("## Done with script ASV_inference.R.\n")
Sys.time()