#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-fq2fa-%j.out

## Report
echo "## Starting script fq2fa.sh"
date
echo

## Software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/seqtk

## Strict settings
set -euo pipefail

## Command-line args
fq_in=$1
fa_out=$2

## Report
echo "## Input FASTQ file:       $fq_in"
echo "## Output FASTA file:      $fa_out"

## Convert FASTQ to FASTA
seqtk seq -a "$fq_in" > "$fa_out"

## Report
echo "## Done with script fq2fa.sh"
date
echo
