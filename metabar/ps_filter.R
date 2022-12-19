#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-ps_filter-%j.out

#TODO - Remove negative controls

# SETUP ------------------------------------------------------------------------
# Packages
packages <- c("BiocManager", "tidyverse", "phyloseq", "decontam",
              "microbiome", "biomformat")

# Default parameter values
contam_method <- "either"     # "prevalence" (neg. control), "frequency" (DNA conc.), "either", "both", or "NA"
contam_thres <- 0.1           # P-value for an ASV to be considered a contaminant 
conc_column <- NA             # Name of the column in the metadata containing DNA concentrations (use 'NA' if none)
batch_column <- NA            # Name of the column in the metadata containing batch IDs (use 'NA' if none)
neg_column <- NA              # Name of the column in the metadata indicating neg. control status; specify either `neg_column` or `neg_ids` to identify negative controls (use 'NA' if none)
neg_ids <- NA                 # IDs of samples that are neg. controls; specify either `neg_column` or `neg_ids` to identify negative controls (use 'NA' if none) 
min_ASV <- 1000               # Min. total ASV count for a sample; sample will be excluded if it has a lower value 
rm_offtarget <- TRUE          # Whether to remove off-target taxa: Chloroplasts, Mitochondria & Eukaryotes

# Process command-line arguments
if (!require(argparse)) install.packages("argparse", repos = "https://cran.rstudio.com/")
library(argparse)
parser <- ArgumentParser()
parser$add_argument("-i", "--ps_in",
                    type = "character", required = TRUE,
                    help = "Input file (phyloseq object RDS) (REQUIRED)")
parser$add_argument("-o", "--ps_out",
                    type = "character", required = TRUE,
                    help = "Output file (phyloseq object RDS) (REQUIRED)")
parser$add_argument("-q", "--outdir_qc",
                    type = "character", required = TRUE,
                    help = "Output dir for QC files (REQUIRED)")
parser$add_argument("-c", "--config",
                    type = "character", default = NULL,
                    help = "Config file (R script)")
args <- parser$parse_args()

ps_in <- args$ps_in
ps_out <- args$ps_out
outdir_qc <- args$outdir_qc
config_file <- args$config

# Load packages
if (!require(pacman)) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(char = packages)

# Source config script
if (! is.null(config_file)) {
  message("\n# Sourcing config file ", config_file)
  source(config_file)
} else {
  message("\n# NOT using a config file")
}

# Check input
if (!is.na(neg_column) & !is.na(neg_ids)) {
    message("\n# ERROR: neg_column & neg_ids can't both be NA")
    message("# neg_column is", neg_column, "and neg_ids is", neg_ids)
    stop()
}

# Define output files
file_id <- sub("ps_", "", sub(".rds", "", basename(ps_out)))
outdir <- dirname(ps_out)
if (is.na(outdir_qc)) outdir_qc <- file.path(outdir, "qc")

outfile_biom <- file.path(outdir, paste0(file_id, ".biom"))
outfile_contam_df <- file.path(outdir_qc, "contam_df.txt")
outprefix_contamplot <- file.path(outdir_qc, "contam_abund-plot")
outfile_contam_plot_df <- file.path(outdir_qc, "contam_abund-plot_df.txt")
outfile_neg_control <- file.path(outdir_qc, "neg_control_ASVs.txt")

# Process parameters
if (is.na(batch_column)) batch_column <- NULL
if (is.na(conc_column)) conc_column <- NULL
if (is.na(neg_column)) neg_column <- NULL

# Create output dirs if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(outdir_qc)) dir.create(outdir_qc, recursive = TRUE)

# Report
message()
message("# ====================================================================")
message("#               STARTING SCRIPT PS_FILTER.R")
message("# ====================================================================")
Sys.time()
message("# Input phyloseq RDS file:                      ", ps_in)
message("# Output phyloseq RDS file:                     ", ps_out)
message("# Output QC file dir:                           ", outdir_qc)
message()
message("# Contaminant-checking method:                  ", contam_method)
message("# Contaminant-checking threshold:               ", contam_thres)
message()
if (!is.null(conc_column)) message("# Metadata column w/ DNA concentrations:        ", conc_column)
if (!is.null(batch_column)) message("# Metadata column w/ batch IDs:                 ", batch_column)
if (!is.null(neg_column)) message("# Metadata column w/ negative controls:         ", neg_column)
if (!is.na(neg_ids)) message("# IDs of samples that are negative controls:    ", neg_ids)
message()
message("# Minimum nr of ASVs for a sample to be kept:   ", min_ASV)
message()
message("# Removing off target taxa?                     ", rm_offtarget)
message("# ====================================================================")
message()


# FUNCTIONS --------------------------------------------------------------------
# Plot of DNA concentration vs ASV abundance
conc_plot <- function(ASVs, plot_id, df) {

    # Add p-value annotation
    df <- df %>%
        arrange(p.freq) %>%
        mutate(
            ASV_p = paste0(ASV, " (p = ", round(p.freq, 4), ")"),
            ASV_p = fct_inorder(ASV_p)
        )

    p <- df %>%
        filter(ASV %in% ASVs) %>%
        ggplot(aes(x = .data[[conc_column]], y = ASV_freq)) +
        geom_smooth(method = "lm", se = FALSE, size = 0.5, color = "red") +
        geom_point() +  #TODO - COLOR BY TREATMENT OR SIMILAR
        facet_wrap(vars(ASV_p)) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x = "DNA concentration", y = "ASV abundance") +
        theme_bw()

    outfile_p_contam <- paste0(outprefix_contamplot, "_", plot_id, ".png")
    ggsave(outfile_p_contam, p, width = 8, height = 8)
}


# READ INPUT DATA --------------------------------------------------------------
ps_raw <- readRDS(ps_in)


# CHECK FOR CONTAMINANTS -------------------------------------------------------
if (!is.na(contam_method)) {
    message("# Identifying and removing potential contaminant ASVs...")

    if (is.null(neg_column) & !is.na(neg_ids)) {
        # If there is no column for neg. control status, create one
        neg_column <- "neg_control"
        sample_data(ps_raw)$neg_control <- sample_names(ps_raw) %in% neg_ids
    }
    if (!is.null(neg_column) & is.na(neg_ids)) {
        # If there is no vector with IDs of neg. control samples, create one
        neg_ids <- sample_names(ps_raw)[which(sample_data(ps_raw)$neg_control == TRUE)]
    }

    if (contam_method != "prevalence") {
      # If contam_method is not prevalence, check DNA concentrations
      # decontam will not accept NA in the DNA concentrations,
      # so we have to remove NAs if DNA concentrations are taken into account

      # Make sure DNA conc of neg control is 0 and not NA:
      neg_idx <- which(sample_names(ps_raw) %in% neg_ids)
      sample_data(ps_raw)[[conc_column]][neg_idx] <- 0.01

      # For contaminant-checking,
      # create a ps object without samples with NA for the DNA concentration
      nsamples_before <- nsamples(ps_raw)

      conc_NAs <- !is.na(sample_data(ps_raw)[[conc_column]])
      ps_raw <- prune_samples(conc_NAs, ps_raw)

      nsamples_after <- nsamples(ps_raw)
      nsamples_removed <- nsamples_before - nsamples_after
      message("# Nr samples w/ no DNA conc excluded for contaminant-checking: ",
              nsamples_removed)
    }

    # Check which ASVs are present in the negative control
    neg_idx <- which(sample_names(ps_raw) %in% neg_ids)
    neg_table <- otu_table(ps_raw)[neg_idx, ] # Subset to neg_control samples
    neg_df <- data.frame(t(neg_table[, which(colSums(neg_table) > 0)])) # Subset to ASVs with >0 counts for neg_control samples
    neg_df <- cbind(neg_df, total_count = rowSums(neg_df)) %>%
        rownames_to_column("ASV") %>%
        arrange(desc(total_count))
    
    message("# Nr distinct ASVs in neg. control(s):  ", nrow(neg_df))
    message("# Total ASV count in neg. control(s):   ", sum(neg_df$total_count))
    write_tsv(neg_df, outfile_neg_control)

    # Check which ASVs are contaminants
    message("# Now running the isContaminant function...")
    contam_df <- isContaminant(ps_raw,
                               method = contam_method,
                               conc = conc_column,
                               neg = neg_column,
                               batch = batch_column,
                               threshold = contam_thres) %>%
        rownames_to_column("ASV")
    write_tsv(contam_df, outfile_contam_df)

    # Report and save basic stats
    n_contam <- sum(contam_df$contaminant)
    message("# FINAL nr of ASVs IDed as contaminants: ", n_contam)

    if (any(contam_df$p.freq < contam_thres, na.rm = TRUE)) {
        # For ASVs IDed as contaminants using the DNA concentration method,
        # create a plot correlating DNA concentration with abundance

        # First get ASVs with a significant value for DNA-conc based testing
        contam_freq_df <- contam_df %>%
            filter(p.freq < contam_thres) %>%
            mutate(ASV = fct_inorder(ASV))

        message("# Nr ASVs IDed as contaminants using DNA conc. only: ",
                nrow(contam_freq_df))

        # Create a df for plotting that includes the ASV abundance values
        contam_plot_df <- otu_table(transform(ps_raw, "compositional")) %>%
            as.data.frame() %>%
            rownames_to_column("sample_id") %>%
            pivot_longer(
                cols = -sample_id,
                names_to = "ASV",
                values_to = "ASV_freq"
            ) %>%
            filter(ASV %in% contam_freq_df$ASV) %>%
            merge(., as(sample_data(ps_raw), "data.frame"),
                  by.x = "sample_id", by.y = "row.names") %>%
            merge(., select(contam_freq_df, ASV, p.freq), by = "ASV") %>%
            mutate(ASV = factor(ASV, levels = levels(contam_freq_df$ASV)))
        
        write_tsv(contam_plot_df, outfile_contam_plot_df)

        # Split the ASVs into groups of 20, and create a plot for each group
        asv_list <- split(contam_freq_df$ASV,
                          ceiling(seq_along(contam_freq_df$ASV) / 12))
        plot_ids <- seq_along(asv_list)
        foo <- mapply(conc_plot, asv_list, plot_ids,
                      MoreArgs = list(df = contam_plot_df))
    }

    if (any(contam_df$p.prev < contam_thres, na.rm = TRUE)) {
        # First get ASVs with a significant value for neg-control based testing
        contam_prev_df <- contam_df %>% filter(p.prev < contam_thres)

        # Check ASVs IDed as contaminants using the negative control
        message("# Nr ASVs IDed as contaminants using neg. control(s) only: ",
                nrow(contam_prev_df))
    }

    if (n_contam > 0) {
        # Remove contaminants
        ps_noncontam <- prune_taxa(!contam_df$contaminant, ps_raw)

        # Process df with info on contaminants
        contam_df$abund <- colSums(ps_raw@otu_table)
        contam_tax <- tax_table(prune_taxa(contam_df$contaminant, ps_raw))
        contam_df_yes <- contam_df %>% filter(contaminant == TRUE)
        contam_df <- merge(contam_df_yes, contam_tax,
                           by.x = "ASV", by.y = "row.names") %>%
            as_tibble() %>%
            arrange(p.freq, p.prev)

        # What proportion of our count data were removed as contaminants?
        prop_kept <- sum(sample_sums(ps_noncontam)) / sum(sample_sums(ps_raw))
        message("# Prop of data kept after contaminant removal: ", prop_kept)
    } else {
        ps_noncontam <- ps_raw
    }
} else {
    ps_noncontam <- ps_raw
}


# REMOVE OFF-TARGET TAXA -------------------------------------------------------
if (rm_offtarget == TRUE) {
  message("\n----------------------------")
  message("# Removing off-target taxa: Chloroplasts, Mitochondria & Eukaryotes...")

  # Create ps subsets for chloroplast, mitochondria, and eukaryotes
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

  # Subset phyloseq object
  bad_taxa <- c(taxa_names(chlr), taxa_names(mit), taxa_names(euk))
  all_taxa <- taxa_names(ps_noncontam)
  good_taxa <- all_taxa[!(all_taxa %in% bad_taxa)]
  ps_target <- prune_taxa(good_taxa, ps_noncontam)

  # What proportion of ASVs were kept?
  prop_kept <- sum(sample_sums(ps_target)) / sum(sample_sums(ps_noncontam))
  n_removed <- nsamples(ps_target) - nsamples(ps_noncontam)
  message("# Prop of data kept after removing off-target taxa:  ", prop_kept)
  message("# Nr of ASVs removed:                                ", n_removed)

  } else {
  message("\n# Not removing off-target taxa...")
  ps_target <- ps_noncontam
}


# FILTER SAMPLES ---------------------------------------------------------------
message("\n------------------------------------")
message("# Removing samples with a total ASV count less than ", min_ASV)
message("# Samples with the lowest total ASV counts:")
sums <- sample_sums(ps_target)
head(sums[order(sums)])

# Remove samples with low counts
ps <- subset_samples(ps_target, sample_sums(ps_target) > min_ASV)
ps <- subset_taxa(ps, taxa_sums(ps) > 0)

# Report how many samples were removed
nsamples_rm <- nsamples(ps_target) - nsamples(ps)
message("# Nr of samples removed after ASV count filtering: ", nsamples_rm)
message("# IDs of samples removed after ASV count filtering:")
setdiff(sample_names(ps_target), sample_names(ps))

# Report how many 0-counts ASVs were removed as consequence of sample filtering
ntaxa_rm <- ntaxa(ps_target) - ntaxa(ps)
message("# Nr of ASVs removed after ASV count filtering: ", ntaxa_rm)
message("# IDs of ASVs removed after ASV count filtering:")
setdiff(taxa_names(ps_target), taxa_names(ps))

# Save RDS file of phyloseq object
saveRDS(ps, ps_out)


# CREATE BIOM FILE -------------------------------------------------------------
ps_biom <- subset_taxa(ps, taxa_sums(ps) > 1) # Filter out singleton ASVs
biom <- make_biom(data = t(otu_table(ps_biom))) # Create biom object
write_biom(biom, outfile_biom)


# WRAP UP ----------------------------------------------------------------------
message("\n---------------------------------")
message("# Number of samples in final object:   ", nsamples(ps))
message("# Number of ASVs in final object:      ", ntaxa(ps))

message("\n# Listing output files:")
system(paste("ls -lh", outfile_contam_df))
if (!is.null(conc_column)) system(paste0("ls -lh ", outprefix_contamplot, "*"))
system(paste("ls -lh", ps_out))

message("\n# Done with script ps_filter.R")
Sys.time()
message()
sessionInfo()
message()
