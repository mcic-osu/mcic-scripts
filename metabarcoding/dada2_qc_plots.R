#!/usr/bin/env Rscript


# SET-UP -----------------------------------------------------------------------

## Process command-line args
args <- commandArgs(trailingOnly = TRUE)
qc_file <- args[1]
outdir <- args[2]

# qc_file <- "results/ASV/main/qc/nseq_summary.txt"
# outdir <- "results/ASV/main/qc"

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

## Load packages
if (! "tidyverse" %in% installed.packages()) install.packages("tidyverse")
library(tidyverse)

## Report
cat("Input file:", qc_file, "\n")
cat("Output dir:", outdir, "\n\n")

## Create output dir if necessary
if (! dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Read QC table
qc <- read.table(qc_file, row.names = 1) %>%
  rownames_to_column("SampleID") %>%
  rename(fastq_filtered = filtered,
         reads_merged = merged,
         non_chimeric = nonchim,
         denoised = denoised_R,
         length_filtered = lenfilter) %>%
  select(-denoised_F)


# PLOTS WITH ABS NUMBERS -------------------------------------------------------

## Barplot
p_bars <- qc %>%
  mutate(fastq_filtering = input - fastq_filtered,
         denoising = fastq_filtered - denoised,
         read_merging = denoised - reads_merged,
         chimera_removal = reads_merged - non_chimeric,
         length_filtering = non_chimeric - length_filtered,
         `(remaining)` = non_chimeric) %>%
  select(SampleID, fastq_filtering, denoising, read_merging,
         chimera_removal, length_filtering, `(remaining)`) %>%
  pivot_longer(cols = -SampleID,
               names_to = "status",
               values_to = "proportion") %>%
  mutate(status = factor(status, levels = status_levels2)) %>%
  ggplot(aes(x = proportion, y = SampleID, fill = status)) +
  geom_col(color = "grey50") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  labs(fill = "Removed by", x = "Proportion of sequences", y = NULL)

## Line plot
p_lines <- qc %>%
  pivot_longer(cols = -SampleID, names_to = "status", values_to = "count") %>%
  mutate(status = factor(status, levels = status_levels)) %>%
  ggplot(aes(y = count, group = SampleID, x = status)) +
  geom_point(size = 1) +
  geom_line(alpha = 0.5) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma) +
  theme_bw() +
  theme(plot.margin = margin(0.2, 1.5, 0.2, 0.2, "cm")) +
  labs(y = "Number of sequences retained")


# PROPORTIONS ------------------------------------------------------------------

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
  select(SampleID, fastq_filtering, denoising, read_merging, chimera_removal,
         length_filtering, `(remaining)`) %>%
  pivot_longer(cols = -SampleID,
               names_to = "status",
               values_to = "proportion") %>%
  mutate(status = factor(status, levels = status_levels2)) %>%
  ggplot(aes(x = proportion, y = SampleID, fill = status)) +
  geom_col(color = "grey50") +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  labs(fill = "Removed by", x = "Proportion of sequences", y = NULL)

p_lines_prop <- qc_prop %>%
  pivot_longer(cols = -SampleID,
               names_to = "status",
               values_to = "proportion_retained") %>%
  mutate(status = factor(status, levels = status_levels)) %>%
  ggplot(aes(y = proportion_retained, group = SampleID, x = status)) +
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
cat("Done with script dada2_qc_plots.R\n")


# POINTS PLOT ------------------------------------------------------------------

# p_points <- qc %>%
#   pivot_longer(cols = -SampleID, names_to = "status", values_to = "count") %>%
#   mutate(status = factor(status, levels = status_levels)) %>%
#   ggplot(aes(x = count, y = SampleID, color = status)) +
#   geom_point(size = 3) +
#   theme_bw() +
#   labs(x = "Number of sequences", y = NULL)

# p_points_prop <- qc_prop %>%
#   pivot_longer(cols = -SampleID,
#                names_to = "status",
#                values_to = "proportion_retained") %>%
#   mutate(status = factor(status, levels = status_levels)) %>%
#   ggplot(aes(x = proportion_retained, y = SampleID, color = status)) +
#   geom_point(size = 3) +
#   scale_x_continuous(limits = c(NA, 1), expand = c(0, 0)) +
#   theme_bw() +
#   labs(y = NULL)

#ggsave(file.path(outdir, "nseq_points.png"), p_points,
#       bg = "white", width = 7, height = 15)

#ggsave(file.path(outdir,"nseq_prop_points.png"), p_points_prop,
#       bg = "white", width = 7, height = 15)
