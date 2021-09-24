#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --output=slurm-fastqc-%j.out

PATH="$HOME/miniconda3/bin:$PATH"
source activate fastqc-env

set -euo pipefail # Strict settings

echo "## Starting script fastqc.sh"
date

## Command-line args:
if [ "$#" -ne 2 ]; then
    echo "Usage: fastqc.sh <input-file.fastq[.gz]> <output-dir>"
    echo "Exiting."
    exit 1
fi

input="$1"
outdir="$2"

## Report:
echo "## Input file:         $input"
echo "## Output directory:   $outdir"

mkdir -p "$outdir"

echo "## Running fastqc..."
fastqc --outdir="$outdir" "$input"

echo -e "\n Done with script."
date
