#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-ps_filter-%j.out

# SETUP ------------------------------------------------------------------------
## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("BiocManager", "tidyverse", "phyloseq", "decontam", "microbiome")
pacman::p_load(char = packages)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ps_in <- args[1]
ps_out <- args[2]
config_file <- args[3]

## Default parameter values
contam_check_method <- "either"  # "prevalence" (neg. control), "frequence" (DNA conc.), "either", "both", or "NA"
contam_check_threshold <- 0.1    # P-value for an ASV to be considered a contaminant
conc_column <- NA                # Name of the column in the metadata containing DNA concentrations
batch_column <- NA               # Name of the column in the metadata containing batch IDs
neg_column <- NA                 # Name of the column in the metadata indicating neg. control status; specify either `neg_column` or `neg_ids` to identify negative controls
neg_ids <- NA                    # IDs of samples that are neg. controls; specify either `neg_column` or `neg_ids` to identify negative controls
min_ASV <- 1000                  # Min. total ASV count for a sample; sample will be excluded if it has a lower value
qc_dir <- NA

## Source config script
source(config_file)

## Define output files
outdir <- dirname(ps_out)
if (is.na(qc_dir)) qc_dir <- file.path(outdir, "qc")
outfile_contam_df <- file.path(qc_dir, "contam_df.txt")
outfile_contam_plot_prefix <- file.path(qc_dir, "contam_abund-plot")

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
Sys.time()
cat("## Input phyloseq RDS file:                        ", ps_in, "\n")
cat("## Output phyloseq RDS file:                       ", ps_out, "\n")
cat("## Output QC file dir:                             ", qc_dir, "\n\n")
cat("## Contaminant-checking method:                    ", contam_check_method, "\n")
cat("## Contaminant-checking threshold:                 ", contam_check_threshold, "\n\n")
cat("## Name of metadata column w/ DNA concentrations:  ", conc_column, "\n")
cat("## Name of metadata column w/ batch IDs:           ", batch_column, "\n")
cat("## Name of metadata column w/ negative controls:   ", neg_column, "\n")
cat("## IDs of samples that are negative controls:      ", neg_ids, "\n\n")
cat("## Minimum nr of ASVs for a sample to be kept:     ", min_ASV, "\n")
cat("-----------------------------\n\n")


# FUNCTIONS --------------------------------------------------------------------
## Plot of DNA concentration vs ASV abundance
conc_plot <- function(ASVs, plot_id, df) {
  
  ## Add p-value annotation
  df <- df %>%
    arrange(p.freq) %>% 
    mutate(ASV_p = paste0(ASV, " (p = ", round(p.freq, 4), ")"),
           ASV_p = fct_inorder(ASV_p))
  
  p <- df %>%
    filter(ASV %in% ASVs) %>%
    ggplot(aes(x = .data[[conc_column]], y = ASV_freq)) +
    geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "red") +
    geom_point(aes(color = salt_level, shape = inoculated)) +
    facet_wrap(vars(ASV_p)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "DNA concentration", y = "ASV abundance") +
    theme_bw()
  
  outfile_contam_plot <- paste0(outfile_contam_plot_prefix, "_", plot_id, ".png")
  ggsave(outfile_contam_plot, p, width = 8, height = 8)
}

# READ INPUT DATA --------------------------------------------------------------
ps_raw <- readRDS(ps_in)


# CHECK FOR CONTAMINANTS -------------------------------------------------------
if ( !is.na(contam_check_method)) {
  cat("## Identifying and removing potential contaminant ASVs...\n")
  
  if (is.null(neg_column) & !is.na(neg_ids)) {
    ## If there is no column for neg. control status, create one
    neg_column <- "neg_control"
    sample_data(ps_raw)$neg_control <- sample_names(ps_raw) %in% neg_ids
  }
  if (!is.null(neg_column) & is.na(neg_ids)) {
    ## If there is no vector with IDs of neg. control samples, create one
    neg_ids <- sample_names(ps_raw)[which(sample_data(ps_raw)$neg_control == TRUE)]
  }
  
  if (contam_check_method != "prevalence") {
    ## decontam will not accept NA in the DNA concentrations,
    ## so we have to remove NAs if DNA concentrations are taken into account
    
    ## Make sure DNA conc of neg control is 0 and not NA:
    neg_idx <- which(sample_names(ps_raw) %in% neg_ids)
    sample_data(ps_raw)[[conc_column]][neg_idx] <- 0.01
    
    ## For contaminant-checking, create a phyloseq object without samples with NA for the DNA concentration
    nsamples_before <- nsamples(ps_raw)
    
    conc_NAs <- !is.na(sample_data(ps_raw)[[conc_column]])
    ps_raw <- prune_samples(conc_NAs, ps_raw)
    
    nsamples_after <- nsamples(ps_raw)
    nsamples_removed <- nsamples_before - nsamples_after
    cat("## For contaminant-checking, excluded", nsamples_removed,
        "samples with no DNA concentration.\n")
  }
  
  ## Check which ASVs are present in the negative control
  # neg_idx <- which(sample_names(ps_raw) %in% neg_ids)
  # neg_table <- otu_table(ps_raw)[neg_idx, ]
  # neg_table[, which(neg_table > 0)]
  
  ## Check which ASVs are contaminants
  cat("## Now running the isContaminant function...\n")
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
      as_tibble() %>%
      arrange(p.freq, p.prev)
    
    cat("## Df with info on ASVs IDed as contaminants:\n")
    print(contam_df)
    write.table(contam_df, outfile_contam_df,
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
    
    if (any(contam_df$p.freq < contam_check_threshold)) {
      ## For ASVs IDed as contaminants using the DNA concentration method,
      ## create a plot correlating DNA concentration with abundance
      
      ## First get the ASVs with a significant value for DNA-conc based testing
      contam_freq_df <- contam_df %>%
        filter(p.freq < contam_check_threshold) %>%
        mutate(ASV = fct_inorder(ASV))
      
      ## Create a df for plotting that includes the ASV aundance values
      contam_plot_df <- otu_table(transform(ps_raw, "compositional")) %>%
        as.data.frame() %>% 
        rownames_to_column("sample_id") %>%
        pivot_longer(cols = -sample_id, names_to = "ASV", values_to = "ASV_freq") %>% 
        filter(ASV %in% contam_freq_df$ASV) %>%
        merge(., as(sample_data(ps_raw), "data.frame"), by = "sample_id") %>%
        merge(., select(contam_freq_df, ASV, p.freq), by = "ASV") %>%
        mutate(ASV = factor(ASV, levels = levels(contam_freq_df$ASV)))
      
      ## Split the ASVs into groups of 20, and create a plot for each group
      ASV_list <- split(contam_freq_df$ASV, ceiling(seq_along(contam_freq_df$ASV) / 12))
      plot_ids <- seq_along(ASV_list)
      foo <- mapply(conc_plot, ASV_list, plot_ids,
                    MoreArgs = list(df = contam_plot_df))
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
cat("\n------------------------------------------\n")
cat("## Removing off-target taxa: Chloroplasts, Mitochondria, and Eukaryotes...\n")

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
cat("\n------------------------------------------\n")
cat("## Removing samples with a total ASV count less than", min_ASV, "...\n")
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
cat("\n------------------------------------------\n")
cat("## Listing output files:\n")
system(paste("ls -lh", ps_out))
system(paste("ls -lh", outfile_contam_df))
if (!is.null(conc_column)) system(paste0("ls -lh ", outfile_contam_plot_prefix, "*"))

cat("\n## Done with script ps_filter.R\n")
Sys.time()
