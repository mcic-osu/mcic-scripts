#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=rcorrector
#SBATCH --output=slurm-rcorr-%j.out

## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate rcorrector-env

## Bash strict mode
set -euo pipefail

## Command-line args
indir="$1"
outdir="$2"
start_at_step="${3-0}"

## Process args
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments - you provided $#" && exit 1
R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')

## Report
echo "## Starting script rcorrector.sh"
date
echo "## Input dir: $indir"
echo "## Output dir: $outdir"
echo
echo "## List of R1 files: $R1_list"
echo "## List of R2 files: $R2_list"
echo
echo "## Start at step: $start_at_step"
echo -e "-------------------------------\n"

## Create output dirs if needed
mkdir -p "$outdir"


# RUN RCORRECTOR ---------------------------------------------------------------
run_rcorrector.pl \
    -t "$SLURM_CPUS_PER_TASK" \
    -od "$outdir" \
    -stage "$start_at_step" \
    -1 "$R1_list" -2 "$R2_list"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing file in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script rcorrector.sh"
date