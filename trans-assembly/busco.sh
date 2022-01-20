#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --job-name=busco
#SBATCH --output=slurm-busco-%j.out


# SETUP ------------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate busco-env

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in=$1
outdir=$2
busco_db=$3

## Check input
[[ "$#" != 3 ]] && echo "## ERROR: Please provide 3 arguments - you provided $#" && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file $fa_in does not exist" && exit 1

## Other paramaters
n_cores=$SLURM_CPUS_PER_TASK

## Report
echo "## Starting script busco.sh"
date
echo "## Input FASTA:         $fa_in"
echo "## Output dir:          $outdir"
echo "## BUSCO db:            $busco_db"
echo "## Number of cores:     $n_cores"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN BUSCO -----------------------------------------------------------------
echo "## Now running busco..."

busco \
    -i "$fa_in" \
    -o "$outdir" \
    -l "$busco_db" \
    -m transcriptome \
    -c "$n_cores" \
    --force


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script busco.sh"
date


# INFO -------------------------------------------------------------------------
#? User guide: https://busco.ezlab.org/busco_userguide.html
#? Conda: https://anaconda.org/bioconda/busco