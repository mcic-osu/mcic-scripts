#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-ps_filter-%j.out

# SETUP ------------------------------------------------------------------------
## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse", "phyloseq", "decontam")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ps_in <- args[1]
ps_out <- args[2]
config_file <- args[3]

#ps_in <- "results/phyloseq/ps_dadatax_raw.rds"
#ps_out <- "results/phyloseq/ps_dadatax_filt.rds"
#config_file <- "workflow/config/ps_filter_config.R"

## Defaults parameters
contam_check_method <- "either"
contam_check_threshold <- 0.1
conc_column <- NA
batch_column <- NA
neg_column <- NA
neg_ids <- NA
min_ASV <- 1000
qc_dir <- NA

## Source config script
source(config_file)

## Define output files
outdir <- dirname(ps_out)
if (is.na(qc_dir)) qc_dir <- file.path(outdir, "qc")
outfile_contam_df <- file.path(qc_dir, "contam_df.txt")
outfile_contam_plot <- file.path(qc_dir, "contam_abund-plot.png")

## Check input
if (!is.na(neg_column) & !is.na(neg_ids)) {
  cat("\n## ERROR: Either neg_column or neg_ids should be something other NA, not both.\n")
  cat("## Currently, neg_column is", neg_column, "and neg_ids is", neg_ids, "\n")
  stop()
}

## Process parameters
if (is.na(batch_column)) batch_column <- NULL
if (is.na(conc_column)) conc_column <- NULL
if (is.na(neg_column)) neg_column <- NULL

## Create output dirs if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

## Report
cat("\n## Starting script ps_filter.R\n")
cat("## Input phyloseq RDS file:", ps_in, "\n")
cat("## Output phyloseq RDS file:", ps_out, "\n")
cat("## Output QC file dir:", qc_dir, "\n\n")
cat("## Contaminant-checking method:", contam_check_method, "\n")
cat("## Contaminant-checking threshold:", contam_check_threshold, "\n\n")
cat("## Name of metadata column containing DNA concentrations:", conc_column, "\n")
cat("## Name of metadata column containing batch IDs:", batch_column, "\n")
cat("## Name of metadata column containing negative controls:", neg_column, "\n")
cat("## IDs of samples that are negative controls:", neg_ids, "\n\n")
cat("## Minimum nr of ASVs for a sample to be kept:", min_ASV, "\n")
cat("-----------------------------\n\n")

# READ INPUT DATA --------------------------------------------------------------
ps_raw <- readRDS(ps_in)


# CHECK FOR CONTAMINANTS -------------------------------------------------------
if ( !is.na(contam_check_method == TRUE)) {
  cat("## Identifying and removing potential contaminant ASVs...\n")
  
  if (is.null(neg_column) & !is.na(neg_ids)) {
    neg_column <- "neg_control"
    sample_data(ps_raw)$neg_control <- sample_names(ps_raw) %in% neg_ids
  }
  if (!is.null(neg_column) & is.na(neg_ids)) {
    neg_ids <- sample_names(ps_raw)[which(sample_data(ps_raw)$neg_control == TRUE)]
  }
  
  if (contam_check_method != "prevalence") {
    ## decontam will not accept NA in the DNA concentrations
    
    ## Make sure DNA conc of neg control is 0:
    neg_idx <- which(sample_names(ps_raw) %in% neg_ids)
    sample_data(ps_raw)[[conc_column]][neg_idx] <- 0.01
    
    ## For contaminant-checking, create a phyloseq object without samples with NA for the DNA concentration
    conc_NAs <- !is.na(sample_data(ps_raw)[[conc_column]])
    ps_raw <- prune_samples(conc_NAs, ps_raw)
  }
  
  ## Check which ASVs are present in the negative control
  # neg_idx <- which(sample_names(ps_raw) %in% neg_ids)
  # neg_table <- otu_table(ps_raw)[neg_idx, ]
  # neg_table[, which(neg_table > 0)]
  
  ## Check which ASVs are contaminants
  contam_df <- isContaminant(ps_raw,
                             method = contam_check_method,
                             conc = conc_column,
                             neg = neg_column,
                             batch = batch_column,
                             threshold = contam_check_threshold)
  
  ## Report and save basic stats
  n_contam <- sum(contam_df$contaminant)
  cat("## Nr of ASVs IDed as contaminants:", n_contam, "\n")
  
  if (n_contam > 0) {
    ## Remove contaminants
    ps_noncontam <- prune_taxa(!contam_df$contaminant, ps_raw)
    
    ## Process df with info on contaminants
    contam_df$abund <- colSums(ps_raw@otu_table)
    contam_tax <- tax_table(prune_taxa(contam_df$contaminant, ps_raw))
    contam_df_yes <- contam_df %>% filter(contaminant == TRUE)
    contam_df <- merge(contam_df_yes, contam_tax, by = "row.names") %>%
      rename(ASV = Row.names) %>%
      as_tibble()
    
    cat("## Stats for ASVs IDed as contaminants:\n")
    print(contam_df)
    write.table(contam_df, outfile_contam_df,
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
    
    if (!is.null(conc_column)) {
      ## Create a plot correlating DNA concentration with abundance for up to the first 6 contaminants
      p_contam <- plot_frequency(ps_raw,
                                 contam_df$ASV,
                                 conc = conc_column) +
        labs(x = "DNA concentration", y = "ASV abundance") +
        theme_bw(base_size = 14)
      ggsave(outfile_contam_plot, p_contam, width = 8, height = 8)
    }
    
    ## What proportion of our count data were removed as contaminants?
    prop_retained <- sum(sample_sums(ps_noncontam)) / sum(sample_sums(ps_raw))
    cat("## Prop count data retained after contaminant removal:", prop_retained, "\n")
    
  } else {
    ps_noncontam <- ps_raw
    }
} else {
  ps_noncontam <- ps_raw
}


# REMOVE OFF-TARGET TAXA -------------------------------------------------------
cat("\n## Removing off-target taxa: Chloroplasts, Mitochondria, and Eukaryotes...\n")

## Create ps subsets for chloroplast, mitochondria, and eukaryotes
if (any(ps_noncontam@tax_table[, "order"] == "Chloroplast", na.rm = TRUE)) {
  chlr <- subset_taxa(ps_noncontam, order == "Chloroplast")
} else {
  chlr <- NULL
}

if (any(ps_noncontam@tax_table[, "family"] == "Mitochondria", na.rm = TRUE)) {
  mit <- subset_taxa(ps_noncontam, family == "Mitochondria")
} else {
  mit <- NULL
}

if (any(ps_noncontam@tax_table[, "domain"] == "Eukaryota", na.rm = TRUE)) {
  euk <- subset_taxa(ps_noncontam, domain == "Eukaryota")
} else {
  euk <- NULL
}

## Subset phyloseq object
bad_taxa <- c(taxa_names(chlr), taxa_names(mit), taxa_names(euk))
all_taxa <- taxa_names(ps_noncontam)
good_taxa <- all_taxa[!(all_taxa %in% bad_taxa)]

ps_target <- prune_taxa(good_taxa, ps_noncontam)

## What proportion of ASVs were kept?
prop_retained <- sum(sample_sums(ps_target)) / sum(sample_sums(ps_noncontam))
cat("## Prop count data retained after removing off-target taxa:", prop_retained, "\n")


# FILTER SAMPLES ---------------------------------------------------------------
cat("\n## Removing samples with a total ASV count less than", min_ASV, "...\n")
cat("## Samples with the lowest total ASV counts:\n")
sums <- sample_sums(ps_target)
head(sums[order(sums)])

## Remove samples with low counts
ps <- subset_samples(ps_target,
                     sample_sums(ps_target) > min_ASV)

## Report how many samples were removed
nsamples_rm <- nrow(ps_target@otu_table) - nrow(ps@otu_table)
cat("## Nr of samples removed after ASV count filtering:", nsamples_rm, "\n")
cat("## IDs of samples removed after ASV count filtering:\n")
setdiff(sample_names(ps_target), sample_names(ps))


# WRAP UP ----------------------------------------------------------------------
## Save RDS file
saveRDS(ps, ps_out)

## Report
cat("\n## Listing output files:\n")
system(paste("ls -lh", ps_out))
system(paste("ls -lh", outfile_contam_df))
if (!is.null(conc_column)) system(paste("ls -lh", outfile_contam_plot))

cat("\n## Done with script ps_filter.R\n")
Sys.time()
