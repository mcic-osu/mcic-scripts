#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --job-name=cdhit
#SBATCH --output=slurm-cdhit-%j.out


# SETUP ------------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate cd-hit-env

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in="$1"
fa_out="$2"

## Process input parameters
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments - you provided $#" && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file $fa_in does not exist" && exit 1

## Other paramaters
n_cores=$SLURM_CPUS_PER_TASK
mem_mb=$((8*(SLURM_MEM_PER_NODE)/10))

## Report
echo "## Starting script cd-hit.sh"
date
echo "## Input FASTA:       $fa_in"
echo "## Output FASTA:      $fa_out"
echo
echo "## Number of cores:   $n_cores"
echo "## Memory in MB:      $mem_mb"
echo -e "-------------------------------\n"

## Create output dirs if needed
outdir=$(dirname "$fa_out")
mkdir -p "$outdir"


# RUN CD-HIT -------------------------------------------------------------------
echo "## Now running CD-HIT-EST..."

cd-hit-est -i "$fa_in" -o "$fa_out" \
    -c 0.95 \
    -n 8 \
    -p 1 \
    -g 1 \
    -M "$mem_mb" \
    -T "$n_cores" \
    -d 40

#? Parameter settings from https://www.protocols.io/view/de-novo-transcriptome-assembly-workflow-ghebt3e?step=10
#? -c similarity threshold
#? -n word_length
#? -p 1 => print alignment overlap
#? -g 1 => slow/accurate cluster mode

# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing the output file:"
ls -lh "$fa_out"

echo -e "\n## Done with script cd-hit.sh"
date
