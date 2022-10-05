#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --output=slurm-pycoqc-%j.out

# SETUP ------------------------------------------------------------------------
## Load software 
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/pycoqc-env

## Bash strict settings
set -ueo pipefail

## Command-line args
infile=$1
outfile=$2

## Input checks
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments; you provided $#" && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input file does not exist" && exit 1

## Create the output directory if it doesn't already exist
outdir=$(dirname "$outfile")
mkdir -p "$outdir"

## Report
echo "## Starting script pycoqc.sh..."
date
echo
echo "## Input FASTQ file:        $infile"
echo "## Output file:             $outfile"
echo -e "--------------------\n"


# RUN PYCOQC -------------------------------------------------------------------
pycoQC \
    -f "$infile" \
    -o "$outfile"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output file:"
ls -lh "$outfile"

echo -e "\n## Done with script pycoqc.sh"
date

## Doc
#? https://tleonardi.github.io/pycoQC/pycoQC/usage/
#? https://tleonardi.github.io/pycoQC/pycoQC/CLI_usage/
