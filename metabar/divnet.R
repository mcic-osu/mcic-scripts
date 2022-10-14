#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-divnet-%j.out

# SET-UP -----------------------------------------------------------------------
## Report
message("\n## Starting script divnet.R")
Sys.time()
message()

## Parse command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i", "--ps_in",
                    type = "character", required = TRUE,
                    help = "Input RDS file with a phyloseq object (REQUIRED)")
parser$add_argument("-o", "--outdir",
                    type = "character", default = "results/divnet",
                    help = "Output dir [default %(default)s]")
parser$add_argument("-f", "--formula",
                    type = "character", default = NULL,
                    help = "R model formula without leading ~ [default %(default)s]")
parser$add_argument("-t", "--tuning",
                    type = "character", default = "fast",
                    help = "Tuning mode ('fast' or 'careful' [default %(default)s]")
parser$add_argument("-b", "--n_boot",
                    type = "integer", default = 2,
                    help = "Number of bootstraps for variance estimation [default %(default)s]")
args <- parser$parse_args()

ps_in <- args$ps_in
outdir <- args$outdir
formula <- args$formula
tuning <- args$tuning
n_boot <- args$n_boot

## Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
if (!require(remotes)) install.packages("remotes", repos = "https://cran.rstudio.com/")
if (!require(breakaway)) remotes::install_github("adw96/breakaway")
if (!require(DivNet)) remotes::install_github("adw96/DivNet")
packages <- c("tidyverse", "phyloseq", "breakaway", "DivNet")
pacman::p_load(char = packages)

## Other parameters
n_cores <- as.integer(system("echo $SLURM_CPUS_PER_TASK", intern = TRUE))

## Report
message("## Input phyloseq RDS file:      ", ps_in)
message("## Output dir:                   ", outdir)
message("## Analysis formula:             ", formula)
message("## MCMC tuning setting:          ", tuning)
message("## Nr of bootstraps:             ", n_boot)
message("## Number of cores:              ", n_cores)
message("------------------------")

## Create output dir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Read input
ps <- readRDS(ps_in)


# PREP PHYLOSEQ OBJECT AND FORMULA ---------------------------------------------
if (!is.null(formula)) {
    ## Infer factors from formula
    factors <- gsub("~", "", formula) %>%
        gsub("\\+", "", .) %>%
        gsub("\\*", "", .) %>%
        gsub(" +", " ", .) %>%
        strsplit(., " ") %>%
        unlist()


    ## Remove samples with NA levels in factors
    message("\n## Nr of samples in original object:     ", nsamples(ps))

    for (factor in factors) {
        message("## Removing NAs for factor:        ", factor)
        ps <- prune_samples(!is.na(sample_data(ps)[[factor]]), ps)
        message("## Nr of samples remaining:        ", nsamples(ps))
    }

    ## Save formula as formula object
    formula <- as.formula(formula)

    file_id <- paste0(factors, collapse = "_")
} else {
    file_id <- "by_ind"
}

    ## Determine the most common ASV #TODO EDIT
    base_tax <- "ASV_1"
    message("## Base taxon:        ", base_tax)

## Define output file
outfile <- file.path(outdir, paste0("divnet_", file_id, ".rds"))


# RUN DIVNET -------------------------------------------------------------------
message("\n## Starting DivNet run...")
#TODO use trycatch first without using taxon base
res <- divnet(ps,
              formula = formula,
              base = base_tax,
              tuning = tuning,
              B = n_boot,
              ncores = n_cores) 


# WRAP UP ----------------------------------------------------------------------
saveRDS(res, outfile)

## Report
message("\n## Listing output file:")
system(paste("ls -lh", outfile))
message("## Done with script divnet.R")
Sys.time()
message()
