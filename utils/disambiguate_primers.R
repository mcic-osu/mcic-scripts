# SETUP ------------------------------------------------------------------------
# Load/install packages
rep <- "https://cloud.r-project.org"
lib <- Sys.getenv("R_LIBS_USER")
dir.create(path = lib, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages( {
  if (!require(argparse)) install.packages("argparse", repos = rep, lib = lib)
  if (!require(pacman)) install.packages("pacman", repos = rep, lib = lib)
  if (!require(BiocManager)) install.packages("BiocManager", repos = rep, lib = lib)
  if (!require(DECIPHER)) BiocManager::install("DECIPHER")
  if (!require(seqinr)) install.packages("seqinr", repos = rep, lib = lib)
  if (!require(purrr)) install.packages("purrr", repos = rep, lib = lib)
} )
packages <- c("argparse", "DECIPHER", "seqinr", "purrr")
pacman::p_load(char = packages, install = TRUE)

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--infile",
                    type = "character", required = TRUE, default = NULL,
                    help = "Input dir with Kallisto output files (REQUIRED)")
parser$add_argument("-o", "--outfile",
                    type = "character", required = FALSE, default = NULL,
                    help = "Output file (default: same as input, with '_disambig' in the name)")
args <- parser$parse_args()
infile <- args$infile
outfile <- args$outfile

# Define/prep outfile
if (!is.null(outfile)) {
  outdir <- dirname(infile)
  outname <- paste0(basename(tools::file_path_sans_ext(infile)),
                    "_disambig.",
                    tools::file_ext(infile))
  outfile <- file.path(outdir, outname)
}
outdir <- dirname(outfile)
dir.create(outdir, showWarnings = FALSE, recursive = FALSE)
stopifnot(outfile != infile)
if (file.exists(outfile)) file.remove(outfile)

# Report
message("\n# Starting script disambiguate_primers.R")
Sys.time()
message()
message("Input file:                     ", infile)
message("Output file:                    ", outfile)
message("======================================================================")
message()

# Functions
write_one <- function(seq_vec, name) {
  message("Primer ", name, " has ", length(seq_vec), " disambiguated sequences")
  if (length(seq_vec) > 1) names <- paste0(name, "_", seq_along(seq_vec))
  if (length(seq_vec) == 1) names <- name 
  walk2(.x = as.character(seq_vec), .y = names, .f = write_fa)
}

write_fa <- function(seq, seqname, filename = outfile) {
  write.fasta(sequences = seq, names = seqname,
              file.out = filename, open = "a",
              as.string = TRUE)
}


# RUN --------------------------------------------------------------------------
# Read the input
seqs <- readDNAStringSet(infile)
disamb_list <- Disambiguate(seqs)

# Run
walk2(.x = disamb_list, .y = names(seqs), .f = write_one)


# WRAP UP ----------------------------------------------------------------------
# List output
message("\n# Listing the output file(s):")
system(paste("ls -lh", outfile))
#system(paste("cat", outfile))
message("\n# Done with script disambiguate_primers.R")
Sys.time()
message()
sessionInfo()
