#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm-metabat-%j.out

## Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate metabat-env

## Bash strict mode
set -euo pipefail

## Command-line args
contigs_in="$1"
bam_in="$2"
prefix_out="$3"

## Other variables
DEPTH_FILE="$prefix_out"_depth.txt

## Create output dir
outdir=$(dirname "$prefix_out")
mkdir -p "$outdir"

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"                            # Retrieve number of cores

## Report
echo "## Starting script metabat.sh..."
date
echo "## Command-line args:"
echo "## Input contigs: $contigs_in"
echo "## Input BAM file: $bam_in"
echo "## Output file prefix: $prefix_out"
echo -e "---------------------------\n\n"

## Get depths
echo "## Running jgi_summarize_bam_contig_depths..."
jgi_summarize_bam_contig_depths \
    --outputDepth "$DEPTH_FILE" \
    "$bam_in"

## Run binning
echo -e "\n## Running metabat..."
metabat \
    -i "$contigs_in" \
    -a "$DEPTH_FILE" \
    -o "$prefix_out" \
    -t "$n_cores" \
    --unbinned

## Report
echo -e "\n## Listing output files:"
ls -lh "$prefix_out"*

echo -e "\n## Done with script metabat.sh"
date