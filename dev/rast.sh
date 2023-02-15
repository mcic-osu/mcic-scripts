#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=rast
#SBATCH --output=slurm-rast-%j.out

# Strict base settings
set -euo pipefail

# Parse command-line args
assembly=$1
outdir=$2
species=$3

# Get the assembly ID
id=$(basename "$assembly" .fasta)

# Report
echo "Starting script rast.sh"
date
echo "Assembly:             $assembly"
echo "Outdir:               $outdir"
echo "Species:              $species"

# Create the output dir
mkdir -pv "$outdir"

# Prepare for annotation by creating a Rast-specific file
echo -e "\n# Creating the RAST 'genome'..."
rast-create-genome \
    --scientific-name "$species" \
    --genetic-code 11 \
    --domain Bacteria \
    --contigs "$assembly" \
    > "$outdir"/"$id".gto

# Run the annotation process
echo -e "\n# Annotating the genome..."
rast-process-genome < "$outdir"/"$id".gto > "$outdir"/"$id"_processed.gto

# Export to several formats
echo -e "\n# Exporting the output..."
rast-export-genome genbank < "$outdir"/"$id"_processed.gto > "$outdir"/"$id".gbk
rast-export-genome gff < "$outdir"/"$id"_processed.gto > "$outdir"/"$id".gff
rast-export-genome feature_data < "$outdir"/"$id"_processed.gto > "$outdir"/"$id".tsv
rast-export-genome contig_fasta < "$outdir"/"$id"_processed.gto > "$outdir"/"$id".fa
rast-export-genome feature_dna < "$outdir"/"$id"_processed.gto > "$outdir"/"$id"_features.fa

# Report
echo "# Listing output files:"
ls -lhd "$PWD"/"$outdir"/*
echo
echo "Done with script"
date

#? Installing the RastTk CLI:
#  - https://www.bv-brc.org/docs/cli_tutorial/cli_installation.html#cli-installation
#  - https://github.com/BV-BRC/BV-BRC-CLI/releases
#  - (See also https://github.com/PATRIC3/PATRIC-distribution/releases/tag/1.039)
# Docs:
#  - https://www.bv-brc.org/docs/cli_tutorial/rasttk_getting_started.html
