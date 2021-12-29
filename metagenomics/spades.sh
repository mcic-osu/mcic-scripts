#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --output=slurm-spades-%j.out

## Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate spades-env

## Bash strict mode
set -euo pipefail

## Command-line args
R1_in="$1"
R2_in="$2"
outdir="$3"

## Create output dir
mkdir -p "$outdir"

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"                            # Retrieve number of cores
mem=$(expr $(expr "$SLURM_MEM_PER_NODE" / 1000) - 1)     # Retrieve memory

## Report
echo "## Starting script spades.sh..."
date
echo "## Command-line args:"
echo "## Input FASTQ file - R1: $R1_in"
echo "## Input FASTQ file - R2: $R2_in"
echo "## Output dir: $outdir"
echo
echo "## Other variables and settings:"
echo "## Number of cores: $n_cores"
echo "## Memory in GB: $mem"
echo -e "---------------------------\n\n"

## Run SPAdes
spades.py \
    --pe1-1 "$R1_in" --pe1-2 "$R2_in" \
    -o "$outdir" \
    -t "$n_cores" -m "$mem" \
    -k 31,51,71,91,111 \
    --meta

## Report
echo "## Done with script spades.sh"
date