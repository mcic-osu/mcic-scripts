#!/usr/bin/env Rscript

# SET-UP -----------------------------------------------------------------------
## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse", "dada2")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
fastq_indir <- args[1]
outdir <- args[2]
config_file <- args[3]            # File (R script) with config for ASV inference
n_cores <- as.integer(args[4])    # Number of computer cores to use

# fastq_indir <- "sandbox/fq_subset"
# outdir <- "results/ASV/subset"
# config_file <- "workflow/config/ASV_config.R"
# n_cores <- 4
# techrep_file <- "metadata/duptable.txt"

## Variable defaults
trunc_f <- 150                    # Truncate F reads after trunc_f bases
trunc_r <- 150                    # Truncate R reads after trunc_r bases
ASV_size_min <- 0                 # Minimum ASV size in bp
ASV_size_max <- Inf               # Minimum ASV size in bp
maxEE <- c(2, 2)                  # Max nr of expected errors in a read
pool <- TRUE                      # Whether or not to using sample pooling in dada algorithm

n_samples <- 5                    # Number of samples to run ("all" => all samples, otherwise specify an integer)
save_rds <- TRUE                  # Whether to save RDS files after every step
techrep_file <- NA
start_at_step <- 1                # At which step to start
                                  # Step 1: Filtering and trimming FASTQ files
                                  # Step 2: Dereplicating FASTQ files
                                  # Step 3: Learning erros
                                  # Step 4: Inferring ASVs
                                  # Step 5: Merging reads
                                  # Step 6: Creating the sequence table
                                  # Step 7: Removing chimeras
                                  # Step 8: Filtering ASVs by length

## Read config file
cat("## Sourcing config file...")
source(config_file)
cat("Done.\n")

## Settings dependent on config
chimera_method <- ifelse(pool == FALSE, "consensus", "pooled")
if (length(maxEE) == 1) maxEE <- c(maxEE, maxEE) # e.g. maxEE of "2" will be turned into "c(2, 2)"

## Define and create output dirs
filter_dir <- file.path(outdir, "fastq_filtered")   # For filtered FASTQ files
rds_dir <- file.path(outdir, "rds_intermed")
qc_dir <- file.path(outdir, "qc")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(filter_dir)) dir.create(filter_dir, recursive = TRUE)
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

## Get paths to input FASTQ files
fq_raw_F <- sort(list.files(fastq_indir, pattern = "_R1_001.fastq.gz",
                            full.names = TRUE))
fq_raw_R <- sort(list.files(fastq_indir, pattern = "_R2_001.fastq.gz",
                            full.names = TRUE))

## If needed, select only a subset of the FASTQ files
n_samples_all <- length(fq_raw_F)
if (n_samples == "all") n_samples <- n_samples_all
fq_raw_F <- fq_raw_F[1:n_samples]
fq_raw_R <- fq_raw_R[1:n_samples]

cat("## Number of samples to be analyzed:", n_samples, "\n")
cat("## First 6 FASTA files:", head(fq_raw_F), "\n")

## Extract sample IDs from FASTQ file names 
sampleIDs <- sub("_L001.*", "", basename(fq_raw_F))

## Define output files
fq_filt_F <- file.path(filter_dir, paste0(sampleIDs, "_F_filt.fastq"))
fq_filt_R <- file.path(filter_dir, paste0(sampleIDs, "_R_filt.fastq"))

errorplot_F_file <- file.path(qc_dir, "errors_F.png")
errorplot_R_file <- file.path(qc_dir, "errors_R.png")

nseq_file <- file.path(qc_dir, "nseq_summary.txt")
nASV_file <- file.path(qc_dir, "nASV_summary.txt")
fasta_out <- file.path(outdir, "ASVs.fa")

seqtab_all_file <- file.path(rds_dir, "seqtab_all.rds")
seqtab_nonchim_file <- file.path(rds_dir, "seqtab_nonchim.rds")
seqtab_lenfilter_file <- file.path(rds_dir, "seqtab_nonchim_lenfilter.rds")
seqtab_final_file <- file.path(outdir, "seqtab.rds")

## Report command-line arguments
cat("## Starting script ASV_inference.R...\n")
Sys.time()
cat("## Dir with input FASTQ files:", fastq_indir, "\n")
cat("## Output dir:", outdir, "\n")
cat("## Truncate F after n bases:", trunc_f, "\n")
cat("## Truncate R after n bases:", trunc_r, "\n")
cat("## Min ASV size:", ASV_size_min, "\n")
cat("## Max ASV size:", ASV_size_max, "\n")
cat("## Max nr of expected errors in a read:", maxEE, "\n")
cat("## Use sample pooling for dada algorithm:", pool, "\n")
cat("## Chimera ID method:", chimera_method, "\n\n")
cat("## Number of cores:", n_cores, "\n")
cat("## Number of samples to analyze:", n_samples, "\n")
cat("## Start at step:", start_at_step, "\n")
cat("## Save RDS files:", save_rds, "\n\n")


# FASTQ FILE FILTERING  --------------------------------------------------------
if (start_at_step <= 1) {
    cat("\n----------------\n## Step 1: Filtering and trimming FASTQ files...\n")

    filter_results <-
      filterAndTrim(fq_raw_F, fq_filt_F,
                    fq_raw_R, fq_filt_R,
                    truncLen = c(trunc_f, trunc_r),
                    trimLeft = 0,
                    trimRight = 0,
                    maxN = 0,
                    maxEE = maxEE,
                    truncQ = 2,
                    rm.phix = TRUE,
                    multithread = n_cores, 
                    compress = FALSE, verbose = TRUE) 

    head(filter_results)

    ## Save objects to RDS files
    if (save_rds) saveRDS(filter_results, file.path(rds_dir, "filter_results.rds"))
}


# FASTQ FILE DEREPLICATION -----------------------------------------------------
if (start_at_step >= 2) {
    filter_results <- readRDS(file.path(rds_dir, "filter_results.rds"))
}

if (start_at_step <= 2) {
    cat("\n----------------\n## Step 2: Dereplicating FASTQ files...\n")

    fq_derep_F <- derepFastq(fq_filt_F, verbose = FALSE)
    fq_derep_R <- derepFastq(fq_filt_R, verbose = FALSE)

    names(fq_derep_F) <- sampleIDs
    names(fq_derep_R) <- sampleIDs

    ## Save objects to RDS files
    if (save_rds) saveRDS(fq_derep_F, file.path(rds_dir, "fq_derep_F.rds"))
    if (save_rds) saveRDS(fq_derep_R, file.path(rds_dir, "fq_derep_R.rds"))
}

# ERROR LEARNING ---------------------------------------------------------------
if (start_at_step >= 3 & start_at_step < 6) {
    fq_derep_F <- readRDS(file.path(rds_dir, "fq_derep_F.rds"))
    fq_derep_R <- readRDS(file.path(rds_dir, "fq_derep_R.rds"))
}

if (start_at_step <= 3) {
    cat("\n----------------\n## Step 3: Learning errors...\n")

    err_F <- learnErrors(fq_derep_F,
                         multithread = n_cores)
    err_R <- learnErrors(fq_derep_R,
                         multithread = n_cores)

    ## Save objects to RDS files
    if (save_rds) saveRDS(err_F, file.path(rds_dir, "err_F.rds"))
    if (save_rds) saveRDS(err_R, file.path(rds_dir, "err_R.rds"))

    ## Plot errors
    p <- plotErrors(err_F, nominalQ = TRUE)
    ggsave(errorplot_F_file, width = 8, height = 7)
    p <- plotErrors(err_R, nominalQ = TRUE)
    ggsave(errorplot_R_file, width = 8, height = 7)
}

# INFER ASVS -------------------------------------------------------------------
if (start_at_step >= 4 & start_at_step < 6) {
    err_F <- readRDS(file.path(rds_dir, "err_F.rds"))
    err_R <- readRDS(file.path(rds_dir, "err_R.rds"))
}

if (start_at_step <= 4) {
    cat("\n----------------\n## Step 4: Inferring ASVs...\n")

    dada_F <- dada(fq_derep_F, err = err_F,
                   pool = pool,
                   multithread = n_cores)
    dada_R <- dada(fq_derep_R, err = err_R,
                   pool = pool,
                   multithread = n_cores)

    ## Save objects to RDS files
    if (save_rds) saveRDS(dada_F, file.path(rds_dir, "dada_F.rds"))
    if (save_rds) saveRDS(dada_R, file.path(rds_dir, "dada_R.rds"))
}

# MERGE READ PAIRS -------------------------------------------------------------
if (start_at_step >= 5) {
    dada_F <- readRDS(file.path(rds_dir, "dada_F.rds"))
    dada_R <- readRDS(file.path(rds_dir, "dada_R.rds"))
}

if (start_at_step <= 5) {
    cat("\n----------------\n## Step 5: Merging read pairs...\n")

    mergers <- mergePairs(dada_F, fq_derep_F,
                          dada_R, fq_derep_R)

    ## Save objects to RDS files
    if (save_rds) saveRDS(mergers, file.path(rds_dir, "mergers.rds"))

    ## Remove large objects from the environment to save memory
    rm(fq_derep_F)
    rm(fq_derep_R)
    foo <- gc()
}


# CREATE SEQUENCE TABLE --------------------------------------------------------
if (start_at_step >= 6) mergers <- readRDS(file.path(rds_dir, "mergers.rds"))

if (start_at_step <= 6) {
    cat("\n----------------\n## Step 6: Creating the sequence table...\n")

    seqtab_all <- makeSequenceTable(mergers)

    ## The dimensions of the object are the nr of samples (rows) and the nr of ASVs (columns):
    cat("## Nr of ASVs before chimera removal:", ncol(seqtab_all), "\n")
    cat("## Total ASV count before chimera removal:", sum(seqtab_all), "\n")

    ## Save objects to RDS files
    if (save_rds) saveRDS(seqtab_all, seqtab_all_file)
}

# REMOVE CHIMERAS --------------------------------------------------------------
if (start_at_step >= 7) seqtab_all <- readRDS(seqtab_all_file)

if (start_at_step <= 7) {
    cat("\n----------------\n## Step 7: Removing chimeras...\n")

    seqtab_nonchim <- removeBimeraDenovo(seqtab_all,
                                         method = chimera_method,
                                         multithread = FALSE) # Multithreading often fails for this step    

    cat("## Nr of ASVs after chimera removal:", ncol(seqtab_nonchim), "\n")
    cat("## Total ASV count after chimera removal:", sum(seqtab_nonchim), "\n")

    ## Save objects to RDS files
    if (save_rds) saveRDS(seqtab_nonchim, seqtab_nonchim_file)
}


# FILTER ASVs BY SIZE  ---------------------------------------------------------
if (start_at_step >= 8) seqtab_nonchim <- readRDS(seqtab_nonchim_file)

cat("## Table of sequence lengths:\n")
print(table(nchar(getSequences(seqtab_nonchim))))

if (start_at_step <= 8 & !(ASV_size_min == 0 & ASV_size_max == Inf)) {
  cat("\n----------------\n## Step 8: Filtering ASVs by length:\n")
  seqtab_lenfilter <- seqtab_nonchim[, nchar(colnames(seqtab_nonchim)) %in% seq(ASV_size_min, ASV_size_max)]

  cat("## Nr of ASVs after filtering by length:", ncol(seqtab_lenfilter), "\n")
  cat("## Total ASV count after filtering by length:", sum(seqtab_lenfilter), "\n")
  cat("## Table of sequence lengths after filtering by length:\n")
  print(table(nchar(getSequences(seqtab_lenfilter))))

  ## Save objects to RDS files
  if (save_rds) saveRDS(seqtab_lenfilter, seqtab_lenfilter_file)
  
  seqtab_current <- seqtab_lenfilter
} else {
  seqtab_current <- seqtab_nonchim
}


# MERGE MULTIPLE ENTRIES FOR INDIVIDUAL SAMPLES --------------------------------
if (!is.na(techrep_file)) {
  cat("\n----------------\n## Step 9: Merging techreps in sequence table...\n")
  
  ## Read the file with technical replicates
  dup_df <- read.table(techrep_file, col.names = c("seqID", "sampleID"))
  n_sample_rep <- length(unique(dup_df$sampleID))
  n_rep <- nrow(dup_df) - n_sample_rep
  cat("## Nr of samples with technical replicates:", n_sample_rep, "\n")
  cat("## Nr of technical replicates:", n_rep, "\n")

  ## Function to merge technical replicates in a sequence table
  sum_dup <- function(sampleID_dup, seqtab, dup_df) {
    seqIDs_dup <- dup_df$seqID[dup_df$sampleID == sampleID_dup]
    dup_idxs <- which(sampleIDs %in% seqIDs_dup)
    stopifnot(length(seqIDs_dup) == length(dup_idxs))
    dup_summed_row <- t(as.data.frame(colSums(seqtab[dup_idxs, ], na.rm = TRUE)))
    rownames(dup_summed_row) <- sampleID_dup
    return(dup_summed_row)
  }
  
  ## Apply sum_dup function to all samples with multiple techreps
  sampleIDs_dup <- unique(dup_df$sampleID) # ID of sample with multiple techreps
  seqtab_dups <- do.call(rbind,
                         lapply(sampleIDs_dup, sum_dup, seqtab_current, dup_df))
  
  ## Create a seqtab with samples with no techreps
  seqtab_nondups <- seqtab_current[! rownames(seqtab_current) %in% dup_df$seqID, ]
  
  ## Merge the two seqtabs (1 with techreps, 1 without techreps) back together
  seqtab_merged <- rbind(seqtab_dups, seqtab_nondups)
  seqtab_merged <- seqtab_merged[order(rownames(seqtab_merged)), ]
  
  ## Check if the total ASV count is the same in the original vs the "merged" seqtab
  stopifnot(sum(seqtab_current) == sum(seqtab_merged))
  ## Check if the number of seqtab rows matches the number of techreps 
  stopifnot(nrow(seqtab_current) - nrow(seqtab_merged) == n_rep)

  cat("## Nr of samples in seqtab before / after merging:",
      nrow(seqtab_current), " / ", nrow(seqtab_merged), "\n")

  seqtab_current <- seqtab_merged
}


# CREATE QC SUMMARY TABLE: NR OF SEQS ------------------------------------------
cat("\n----------------\n## Creating QC summary tables...\n")

## Define a function to get the nr of unique sequences
getN <- function(x) sum(getUniques(x))

## Calculate nr of unique sequences across denoised and merged seqs
denoised_F <- sapply(dada_F, getN)
denoised_R <- sapply(dada_R, getN)
merged <- sapply(mergers, getN)

## Put together the final QC table
nseq_summary <- data.frame(filter_results,
                           denoised_F,
                           denoised_R,
                           merged,
                           nonchim = rowSums(seqtab_nonchim),
                           lenfilter = rowSums(seqtab_lenfilter),
                           row.names = sampleIDs)
colnames(nseq_summary)[1:2] <- c("input", "filtered")

## Have a look at the first few rows
cat("## First few rows of the QC table with nr of sequences:\n")
head(nseq_summary)

## Write QC summary table to file
write.table(nseq_summary, file = nseq_file,
            sep = "\t", quote = FALSE, row.names = TRUE)


# CREATE QC SUMMARY TABLE: NR OF UNIQUE ASVs -----------------------------------
nASV_summary <- data.frame(step = c("initial", "non_chimeric", "lenfilter"),
                           n_asv = c(ncol(seqtab_all), ncol(seqtab_nonchim),
                                     ncol(seqtab_lenfilter)))

write.table(nASV_summary, file = nASV_file,
            sep = "\t", quote = FALSE, row.names = TRUE)


# CREATE AND WRITE FINAL SEQTAB AND FASTA FILE ---------------------------------
## Save final seqtab RDS
saveRDS(seqtab_current, seqtab_final_file)

## Prepare sequences and headers
asv_seqs <- colnames(seqtab_current)
asv_headers <- paste(">ASV", 1:ncol(seqtab_current), sep = "_")

## Interleave headers and sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))

## Write FASTA file
write(asv_fasta, file = fasta_out)


# LIST OUTPUT FILES ------------------------------------------------------------
cat("\n----------------\n## Listing output files:\n")

cat("## First few filtered FASTQ files:\n", fq_filt_F[1:2], "\n")

cat("## Error profile plot - F:", errorplot_F_file, "\n")
cat("## Error profile plot - R:", errorplot_R_file, "\n")

cat("## FASTA:", fasta_out, "\n")
cat("## Sequence table:", seqtab_final_file, "\n")
cat("## QC table:", nseq_file, "\n")

cat("## Done with script ASV_inference.R.\n")
Sys.time()