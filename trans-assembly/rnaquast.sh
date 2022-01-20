#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --job-name=rnaquast
#SBATCH --output=slurm-rnaquast-%j.out


# SETUP ------------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate rnaquast-env

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in=$1
outdir=$2

## Check input
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments - you provided $#" && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file $fa_in does not exist" && exit 1

## Other paramaters
n_cores=$SLURM_CPUS_PER_TASK

## Report
echo "## Starting script rnaquast.sh"
date
echo "## Input FASTA:         $fa_in"
echo "## Output dir:          $outdir"
echo "## Number of cores:     $n_cores"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN RNAQUAST -----------------------------------------------------------------
echo "## Now running rnaQUAST..."

rnaQUAST.py --transcripts "$fa_in" \
    --output_dir "$outdir" \
    --threads "$n_cores" \
    --strand_specific \
    --gene_mark

#? Not running BUSCO with rnaQUAST because it failed

# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script rnaquast.sh"
date


#? rnaQUAST website: https://github.com/ablab/rnaquast