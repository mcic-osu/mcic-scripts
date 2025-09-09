#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=astral
#SBATCH --output=slurm-astral-%j.out

#? ASTRAL Docs: https://github.com/smirarab/ASTRAL
#? ASTRAL Tutorial: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md

## Command-line args
infile=$1
outfile=$2

## Load software
module load miniconda3/4.12.0-py39
source activate /fs/ess/PAS0471/jelmer/conda/astral-5.7.8

## Strict bash settings
set -euo pipefail

## Report
echo "Starting script astral.sh"
date

## Run ASTRAL
astral -i "$infile" -o "$outfile"

echo "Done with script"
date
