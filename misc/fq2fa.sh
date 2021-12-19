#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-fq2fa-%j.out

## Software
source ~/.bashrc
source activate seqtk-env

## Strict settings
set -euo pipefail

## Command-line args
fq_in=$1
fa_out=$2

## Report
echo "## Starting script trimmomatic.sh"
date
echo "## Input FASTQ file: $fq_in"
echo "## Output FASTA file: $fa_out"

## Convert FASTQ to FASTA
seqtk seq -a "$fq_in" > "$fa_out"

## Report
echo "Done with script."
date