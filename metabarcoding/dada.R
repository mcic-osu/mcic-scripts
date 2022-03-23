#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --output=slurm-dada-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

# SET-UP -----------------------------------------------------------------------
## Report
message("## Starting script dada.R")
Sys.time()
message()

## Parse command-line arguments
if(!"argparse" %in% installed.packages()) install.packages("argparse")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--indir",
                    type = "character", required = TRUE,
                    help = "Input dir with FASTQ files (REQUIRED)")
parser$add_argument("-o", "--outdir",
                    type = "character", default = "results/dada",
                    help = "Output directory [default %(default)s]")
parser$add_argument("-c", "--config",
                    type = "character", default = NULL,
                    help = "Config file")
parser$add_argument("-t", "--threads",
                    type = "integer", default = NULL,
                    help = "Number of threads [default %(default)s]")
args <- parser$parse_args()

fastq_indir <- args$indir
outdir <- args$outdir
config_file <- args$config
n_threads <- args$threads

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("BiocManager", "tidyverse", "dada2")
pacman::p_load(char = packages)

## Other variables
if (is.null(n_threads)) {
    n_threads <- as.integer(system("echo $SLURM_CPUS_PER_TASK", intern = TRUE))
}

## Variable defaults
trunc_f <- 150           # Truncate F reads after trunc_f bases
trunc_r <- 150           # Truncate R reads after trunc_r bases
asv_size_min <- 0        # Minimum ASV size in bp
asv_size_max <- Inf      # Minimum ASV size in bp
max_ee <- c(2, 2)        # Max nr of expected errors in a read
pool <- TRUE             # Whether or not to using sample pooling in dada algorithm

n_samples <- "all"       # Nr of samples to run ("all" => all samples, or specify an integer)
save_rds <- TRUE         # Whether to save RDS files after every step
techrep_file <- NA       # File with technical replicates
start_at_step <- 1       # At which step to start

# Step 1: Filtering and trimming FASTQ files
# Step 2: Dereplicating FASTQ files
# Step 3: Learning erros
# Step 4: Inferring ASVs
# Step 5: Merging reads
# Step 6: Creating the sequence table
# Step 7: Removing chimeras
# Step 8: Filtering ASVs by length

## Read config file
if (! is.null(config_file)) {
  message("## Sourcing config file...")
  source(config_file)
}

## Settings dependent on config
chimera_method <- ifelse(pool == FALSE, "consensus", "pooled")

## Define and create output dirs
filter_dir <- file.path(outdir, "fastq_filtered")     # For filtered FASTQ files
rds_dir <- file.path(outdir, "rds_intermed")          # For intermediate RDS files
qc_dir <- file.path(outdir, "qc")                     # For QC plots and tables

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(filter_dir)) dir.create(filter_dir, recursive = TRUE)
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

## Get paths to input FASTQ files
fq_raw_f <- sort(list.files(
    fastq_indir, pattern = "_R1_001.fastq.gz", full.names = TRUE
))
fq_raw_r <- sort(list.files(
    fastq_indir, pattern = "_R2_001.fastq.gz", full.names = TRUE
))

## If needed, select only a subset of the FASTQ files
n_samples_all <- length(fq_raw_f)
if (n_samples == "all") n_samples <- n_samples_all
fq_raw_f <- fq_raw_f[1:n_samples]
fq_raw_r <- fq_raw_r[1:n_samples]

## Extract sample IDs from FASTQ file names
sample_ids <- sub("(_S\\d\\d)?_L00\\d_R1_.*", "", basename(fq_raw_f))

## Define output files
fq_filt_f <- file.path(filter_dir, paste0(sample_ids, "_f_filt.fastq"))
fq_filt_r <- file.path(filter_dir, paste0(sample_ids, "_r_filt.fastq"))

errorplot_f_file <- file.path(qc_dir, "errors_f.png")
errorplot_r_file <- file.path(qc_dir, "errors_r.png")

nseq_file <- file.path(qc_dir, "nseq_summary.txt")
nasv_file <- file.path(qc_dir, "nasv_summary.txt")
fasta_out <- file.path(outdir, "ASVs.fa")

seqtab_all_file <- file.path(rds_dir, "seqtab_all.rds")
seqtab_nochim_file <- file.path(rds_dir, "seqtab_nochim.rds")
seqtab_lenfilt_file <- file.path(rds_dir, "seqtab_nochim_lenfilter.rds")
seqtab_final_file <- file.path(outdir, "seqtab.rds")

## Report command-line arguments
message("## Dir with input FASTQ files:             ", fastq_indir)
message("## Output dir:                             ", outdir)
message("## Truncate F after n bases:               ", trunc_f)
message("## Truncate R after n bases:               ", trunc_r)
message("## Min ASV size:                           ", asv_size_min)
message("## Max ASV size:                           ", asv_size_max)
cat("## Max nr expected errors per read:            ", max_ee, "\n")
message("## Use sample pooling for dada algorithm:  ", pool)
message("## Chimera ID method:                      ", chimera_method)
message()
message("## Number of cores:                        ", n_threads)
message("## Start at step:                          ", start_at_step)
message("## Save RDS files:                         ", save_rds)
message()
message("## Number of samples to analyze:           ", n_samples)
cat("## First 2 FASTQ files:", head(fq_raw_f, 2), "\n")
cat("## First 2 sample IDs:", head(sample_ids, 2), "\n")
message()


# FASTQ FILE FILTERING  --------------------------------------------------------
if (start_at_step <= 1) {
    message("\n--------------\n## Step 1: Filtering and trimming FASTQ files...")

    filt_res <- filterAndTrim(fq_raw_f, fq_filt_f,
                              fq_raw_r, fq_filt_r,
                              truncLen = c(trunc_f, trunc_r),
                              trimLeft = 0, trimRight = 0,
                              maxN = 0, maxEE = max_ee,
                              truncQ = 2,
                              rm.phix = TRUE,
                              multithread = n_threads,
                              compress = FALSE, verbose = TRUE)
    head(filt_res)

    ## Save objects to RDS files
    if (save_rds) saveRDS(filt_res, file.path(rds_dir, "filt_res.rds"))
}


# FASTQ FILE DEREPLICATION -----------------------------------------------------
if (start_at_step >= 2) {
    filt_res <- readRDS(file.path(rds_dir, "filt_res.rds"))
}

if (start_at_step <= 2) {
    message("\n----------------\n## Step 2: Dereplicating FASTQ files...")

    fq_derep_f <- derepFastq(fq_filt_f, verbose = FALSE)
    fq_derep_r <- derepFastq(fq_filt_r, verbose = FALSE)

    names(fq_derep_f) <- sample_ids
    names(fq_derep_r) <- sample_ids

    ## Save objects to RDS files
    if (save_rds) saveRDS(fq_derep_f, file.path(rds_dir, "fq_derep_f.rds"))
    if (save_rds) saveRDS(fq_derep_r, file.path(rds_dir, "fq_derep_r.rds"))
}

# ERROR LEARNING ---------------------------------------------------------------
if (start_at_step >= 3 & start_at_step < 6) {
    fq_derep_f <- readRDS(file.path(rds_dir, "fq_derep_f.rds"))
    fq_derep_r <- readRDS(file.path(rds_dir, "fq_derep_r.rds"))
}

if (start_at_step <= 3) {
    message("\n----------------\n## Step 3: Learning errors...")

    err_f <- learnErrors(fq_derep_f, multithread = n_threads)
    err_r <- learnErrors(fq_derep_r, multithread = n_threads)

    ## Save objects to RDS files
    if (save_rds) saveRDS(err_f, file.path(rds_dir, "err_f.rds"))
    if (save_rds) saveRDS(err_r, file.path(rds_dir, "err_r.rds"))

    ## Plot errors
    p <- plotErrors(err_f, nominalQ = TRUE)
    ggsave(errorplot_f_file, width = 8, height = 7)
    p <- plotErrors(err_r, nominalQ = TRUE)
    ggsave(errorplot_r_file, width = 8, height = 7)
}

# INFER ASVS -------------------------------------------------------------------
if (start_at_step >= 4 & start_at_step < 6) {
    err_f <- readRDS(file.path(rds_dir, "err_f.rds"))
    err_r <- readRDS(file.path(rds_dir, "err_r.rds"))
}

if (start_at_step <= 4) {
    message("\n----------------\n## Step 4: Inferring ASVs...")

    dada_f <- dada(fq_derep_f, err = err_f, pool = pool,
                   multithread = n_threads)
    dada_r <- dada(fq_derep_r, err = err_r, pool = pool,
                   multithread = n_threads)

    ## Save objects to RDS files
    if (save_rds) saveRDS(dada_f, file.path(rds_dir, "dada_f.rds"))
    if (save_rds) saveRDS(dada_r, file.path(rds_dir, "dada_r.rds"))
}

# MERGE READ PAIRS -------------------------------------------------------------
if (start_at_step >= 5) {
    dada_f <- readRDS(file.path(rds_dir, "dada_f.rds"))
    dada_r <- readRDS(file.path(rds_dir, "dada_r.rds"))
}

if (start_at_step <= 5) {
    message("\n----------------\n## Step 5: Merging read pairs...")

    mergers <- mergePairs(dada_f, fq_derep_f,
                          dada_r, fq_derep_r)

    ## Save objects to RDS files
    if (save_rds) saveRDS(mergers, file.path(rds_dir, "mergers.rds"))

    ## Remove large objects from the environment to save memory
    rm(fq_derep_f)
    rm(fq_derep_r)
    foo <- gc()
}


# CREATE SEQUENCE TABLE --------------------------------------------------------
if (start_at_step >= 6) mergers <- readRDS(file.path(rds_dir, "mergers.rds"))

if (start_at_step <= 6) {
    message("\n----------------\n## Step 6: Creating the sequence table...")

    seqtab_all <- makeSequenceTable(mergers)

    message("## Nr of ASVs before chimera removal:      ", ncol(seqtab_all))
    message("## Total ASV count before chimera removal: ", sum(seqtab_all))

    ## Save objects to RDS files
    if (save_rds) saveRDS(seqtab_all, seqtab_all_file)
}

# REMOVE CHIMERAS --------------------------------------------------------------
if (start_at_step >= 7) seqtab_all <- readRDS(seqtab_all_file)

if (start_at_step <= 7) {
    message("\n----------------\n## Step 7: Removing chimeras...")

    seqtab_nochim <- removeBimeraDenovo(seqtab_all,
                                        method = chimera_method,
                                        multithread = FALSE)
    # Multithreading often fails for this step

    message("## Nr of ASVs after chimera removal:       ", ncol(seqtab_nochim))
    message("## Total ASV count after chimera removal:  ", sum(seqtab_nochim))

    ## Save objects to RDS files
    if (save_rds) saveRDS(seqtab_nochim, seqtab_nochim_file)
}


# FILTER ASVs BY SIZE  ---------------------------------------------------------
if (start_at_step >= 8) seqtab_nochim <- readRDS(seqtab_nochim_file)

message("## Table of sequence lengths:")
print(table(nchar(getSequences(seqtab_nochim))))

if (start_at_step <= 8 & !(asv_size_min == 0 & asv_size_max == Inf)) {
    message("\n----------------\n## Step 8: Filtering ASVs by length:")
    
    asv_range <- nchar(colnames(seqtab_nochim))
    asv_tres <- seq(asv_size_min, asv_size_max)
    seqtab_lenfilt <- seqtab_nochim[, asv_range %in% asv_tres]

    message("## Nr of ASVs after length filtering:       ", ncol(seqtab_lenfilt))
    message("## Total ASV count after length filtering:  ", sum(seqtab_lenfilt))
    message("## Table of sequence lengths after filtering by length:")
    print(table(nchar(getSequences(seqtab_lenfilt))))

    ## Save objects to RDS files
    if (save_rds) saveRDS(seqtab_lenfilt, seqtab_lenfilt_file)

    seqtab_curr <- seqtab_lenfilt
} else {
    seqtab_curr <- seqtab_nochim
}


# MERGE MULTIPLE ENTRIES FOR INDIVIDUAL SAMPLES --------------------------------
if (!is.na(techrep_file)) {
    message("\n-------------\n## Step 9: Merging techreps in sequence table...")

    ## Read the file with technical replicates
    dup_df <- read.table(techrep_file, col.names = c("seqID", "sample_id"))
    n_sample_rep <- length(unique(dup_df$sample_id))
    n_rep <- nrow(dup_df) - n_sample_rep
    message("## Nr of samples with technical replicates:  ", n_sample_rep)
    message("## Nr of technical replicates:               ", n_rep)

    ## Function to merge technical replicates in a sequence table
    sum_dup <- function(sampleid_dup, seqtab, dup_df) {
        seqids_dup <- dup_df$seqID[dup_df$sample_id == sampleid_dup]
        dup_idxs <- which(sample_ids %in% seqids_dup)
        stopifnot(length(seqids_dup) == length(dup_idxs))
        dup_sum_row <- t(as.data.frame(colSums(seqtab[dup_idxs, ],
                                               na.rm = TRUE)))
        rownames(dup_sum_row) <- sampleid_dup
        return(dup_sum_row)
    }

    ## Apply sum_dup function to all samples with multiple techreps
    sample_ids_dup <- unique(dup_df$sample_id) # ID of sample with techreps
    seqtab_dups <- do.call(
        rbind,
        lapply(sample_ids_dup, sum_dup, seqtab_curr, dup_df)
    )

    ## Create a seqtab with samples with no techreps
    seqtab_nondups <- seqtab_curr[!rownames(seqtab_curr) %in% dup_df$seqID, ]

    ## Merge the two seqtabs (1 with techreps, 1 without techreps) back together
    seqtab_merged <- rbind(seqtab_dups, seqtab_nondups)
    seqtab_merged <- seqtab_merged[order(rownames(seqtab_merged)), ]

    ## Check if total ASV count is same in the original vs the "merged" seqtab
    stopifnot(sum(seqtab_curr) == sum(seqtab_merged))
    ## Check if the number of seqtab rows matches the number of techreps
    stopifnot(nrow(seqtab_curr) - nrow(seqtab_merged) == n_rep)

    message("## Nr of samples in seqtab before / after merging:",
            nrow(seqtab_curr), " / ", nrow(seqtab_merged))

    seqtab_curr <- seqtab_merged
}


# CREATE QC SUMMARY TABLE: NR OF SEQS ------------------------------------------
message("\n----------------\n## Creating QC summary tables...")

## Calculate nr of unique sequences across denoised and merged seqs
denoised_f <- sapply(dada_f, function(x) sum(getUniques(x)))
denoised_r <- sapply(dada_r, function(x) sum(getUniques(x)))
merged <- sapply(mergers, function(x) sum(getUniques(x)))

## Put together the final QC table
nseq_summary <- data.frame(
    sample_id = sample_ids,
    filt_res,
    denoised_f,
    denoised_r,
    merged,
    nonchim = rowSums(seqtab_nochim),
    row.names = NULL
)
if (!(asv_size_min == 0 & asv_size_max == Inf)) {
    nseq_summary$lenfilter <- rowSums(seqtab_lenfilt)
}
colnames(nseq_summary)[2:3] <- c("input", "filtered")
write_tsv(nseq_summary, nseq_file)

## Have a look at the first few rows
message("## First few rows of the QC table with nr of sequences:")
head(nseq_summary)


# CREATE QC SUMMARY TABLE: NR OF UNIQUE ASVs -----------------------------------
nasv_summary <- data.frame(step = c("initial", "non_chimeric"),
                           n_asv = c(ncol(seqtab_all), ncol(seqtab_nochim)))
if (!(asv_size_min == 0 & asv_size_max == Inf)) {
    nasv_summary$lenfilter <- ncol(seqtab_lenfilt)
}
write_tsv(nasv_summary, nasv_file)


# CREATE AND WRITE FINAL SEQTAB AND FASTA FILE ---------------------------------
## Save final seqtab RDS
saveRDS(seqtab_curr, seqtab_final_file)

## Prepare sequences and headers
asv_seqs <- colnames(seqtab_curr)
asv_headers <- paste(">ASV", seq_len(ncol(seqtab_curr)), sep = "_")

## Interleave headers and sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))

## Write FASTA file
write(asv_fasta, file = fasta_out)


# LIST OUTPUT FILES ------------------------------------------------------------
message("\n-----------------------\n## Listing output files...")
message("## First filtered FASTQ file:")
system(paste("ls -lh", fq_filt_f[1]))
message("## Plot with F reads error profile:")
system(paste("ls -lh", errorplot_f_file))
message("## FASTA file:")
system(paste("ls -lh", fasta_out))
message("## Sequence table:")
system(paste("ls -lh", seqtab_final_file))
message("## QC table -- reads:")
system(paste("ls -lh", nseq_file))
message("## QC table -- ASVs:")
system(paste("ls -lh", nasv_file))

message("\n## Done with script dada.R.")
Sys.time()
message()
