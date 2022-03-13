#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-dada2-qc-plots-%j.out

# SET-UP -----------------------------------------------------------------------
message("## Starting script ASV_qcplots.R")
Sys.time()
message()

## Parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--infile",
                    type = "character", required = TRUE,
                    help = "Input file with dada2 QC stats (REQUIRED)")
parser$add_argument("-o", "--outdir",
                    type = "character", default = "results/dada/qc",
                    help = "Output directory [default %(default)s]")
args <- parser$parse_args()

infile <- args$infile
outdir <- args$outdir

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse")
pacman::p_load(char = packages)

## Report
message("## Input file:     ", infile)
message("## Output dir:     ", outdir)
message()

## Files and settings
props_df_file <- file.path(outdir, "nseq_props.txt")
meanprops_df_file <- file.path(outdir, "nseq_meanprops.txt")
meancounts_df_file <- file.path(outdir, "nseq_mean.txt")

plotfile_bars <- file.path(outdir, "nseq_bars.png")
plotfile_lines <- file.path(outdir,"nseq_lines.png")
plotfile_bars_prop <- file.path(outdir,"nseq_prop_bars.png")
plotfile_lines_prop <- file.path(outdir, "nseq_prop_lines.png")

### Plotting order
status_levels <- c("input", "fastq_filtered", "denoised", "reads_merged",
                   "non_chimeric", "length_filtered")
status_levels2 <- c("fastq_filtering", "denoising", "read_merging",
                    "chimera_removal", "length_filtering", "(remaining)")

## Create output dir if necessary
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# READ INPUT FILES -------------------------------------------------------------
qc <- read_tsv(infile, show_col_types = FALSE) %>%
  rename(fastq_filtered = filtered,
         reads_merged = merged,
         non_chimeric = nonchim,
         denoised = denoised_r,
         length_filtered = lenfilter) %>%
  select(-denoised_f)
colnames(qc)[1] <- "sample_id"

# PLOTS WITH ABSOLUTE NUMBERS --------------------------------------------------
## Barplot
p_bars <- qc %>%
  mutate(fastq_filtering = input - fastq_filtered,
         denoising = fastq_filtered - denoised,
         read_merging = denoised - reads_merged,
         chimera_removal = reads_merged - non_chimeric,
         length_filtering = non_chimeric - length_filtered,
         `(remaining)` = length_filtered) %>%
  select(sample_id, fastq_filtering, denoising, read_merging,
         chimera_removal, length_filtering, `(remaining)`) %>%
  pivot_longer(cols = -sample_id,
               names_to = "status",
               values_to = "proportion") %>%
  mutate(status = factor(status, levels = status_levels2)) %>%
  ggplot(aes(x = proportion, y = sample_id, fill = status)) +
  geom_col(color = "grey50") +
  scale_x_continuous(expand = c(0, 0), labels = scales::comma) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  labs(fill = "Removed by", x = "Number of sequences", y = NULL)

## Line plot
p_lines <- qc %>%
  pivot_longer(cols = -sample_id, names_to = "status", values_to = "count") %>%
  mutate(status = factor(status, levels = status_levels)) %>%
  ggplot(aes(y = count, group = sample_id, x = status)) +
  geom_point(size = 1) +
  geom_line(alpha = 0.5) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma) +
  theme_bw() +
  theme(plot.margin = margin(0.2, 1.5, 0.2, 0.2, "cm")) +
  labs(y = "Number of sequences retained")


# PLOTS WITH PROPORTIONS -------------------------------------------------------
## Make df with proportions
qc_prop <- qc %>%
  mutate(fastq_filtered = round(fastq_filtered / input, 3),
         denoised = round(denoised / input, 3),
         reads_merged = round(reads_merged / input, 3),
         non_chimeric = round(non_chimeric / input, 3),
         length_filtered = round(length_filtered / input, 3),
         input = 1)

## Create plots
p_bars_prop <- qc_prop %>%
  mutate(fastq_filtering = input - fastq_filtered,
         denoising = fastq_filtered - denoised,
         read_merging = denoised - reads_merged,
         chimera_removal = reads_merged - non_chimeric,
         length_filtering = non_chimeric - length_filtered,
         `(remaining)` = length_filtered) %>%
  select(sample_id, fastq_filtering, denoising, read_merging, chimera_removal,
         length_filtering, `(remaining)`) %>%
  pivot_longer(cols = -sample_id,
               names_to = "status",
               values_to = "proportion") %>%
  mutate(status = factor(status, levels = status_levels2)) %>%
  ggplot(aes(x = proportion, y = sample_id, fill = status)) +
  geom_col(color = "grey50") +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  labs(fill = "Removed by", x = "Proportion of sequences", y = NULL)

p_lines_prop <- qc_prop %>%
  pivot_longer(cols = -sample_id,
               names_to = "status",
               values_to = "proportion_retained") %>%
  mutate(status = factor(status, levels = status_levels)) %>%
  ggplot(aes(y = proportion_retained, group = sample_id, x = status)) +
  geom_point(size = 1) +
  geom_line(alpha = 0.5) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(limits = c(NA, 1), expand = c(0, 0)) +
  theme_bw() +
  theme(plot.margin = margin(0.2, 1.5, 0.2, 0.2, "cm")) +
  labs(y = "Proportion of sequences retained")

## Get mean proportions
qc_prop_mean <- qc_prop %>% summarise_if(is.numeric, ~ round(mean(.x), 4))
qc_mean <- qc %>% summarise_if(is.numeric, ~ round(mean(.x)))


# SAVE PLOTS AND DFs -----------------------------------------------------------
ggsave(plotfile_bars, p_bars,
       bg = "white", width = 7, height = 15)
ggsave(plotfile_lines, p_lines,
       bg = "white", width = 7, height = 5)

ggsave(plotfile_bars_prop, p_bars_prop,
       bg = "white", width = 7, height = 15)
ggsave(plotfile_lines_prop, p_lines_prop,
       bg = "white", width = 7, height = 5)

write_tsv(qc_prop, props_df_file)
write_tsv(qc_prop_mean, meanprops_df_file)
write_tsv(qc_mean, meancounts_df_file)

## Report
message("\n## Listing output files:")
system(paste("ls -lh", props_df_file))
system(paste("ls -lh", meanprops_df_file))
system(paste("ls -lh", meancounts_df_file))
system(paste("ls -lh", plotfile_bars))
system(paste("ls -lh", plotfile_lines))
system(paste("ls -lh", plotfile_bars_prop))
system(paste("ls -lh", plotfile_lines_prop))

message("\n## Done with script ASV_qcplots.R")
Sys.time()
message()
