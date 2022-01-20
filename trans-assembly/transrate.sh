#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --job-name=transrate
#SBATCH --output=slurm-transrate-%j.out


# SETUP ------------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate transrate-env

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in=$1
fq_dir=$2
outdir=$3

## Check input
[[ "$#" != 3 ]] && echo "## ERROR: Please provide 3 arguments - you provided $#" && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file $fa_in does not exist" && exit 1
[[ ! -d "$fq_dir" ]] && echo "## ERROR: Input dir $fq_dir does not exist" && exit 1

## Process input parameters
R1_list=$(echo "$fq_dir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$fq_dir"/*R2*fastq.gz | sed 's/ /,/g')

## Other paramaters
n_cores=$SLURM_CPUS_PER_TASK

## Report
echo "## Starting script transrate.sh"
date
echo "## Input FASTA:         $fa_in"
echo "## Input FASTQ:         $fq_dir"
echo "## Output dir:          $outdir"
echo
echo "## Number of cores:     $n_cores"
echo
echo "## List of R1 files:    $R1_list"
echo "## List of R2 files:    $R2_list"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN TRANSRATE ----------------------------------------------------------------
echo "## Now running TransRate..."

transrate --assembly "$fa_in" -o "$outdir" \
    --left "$R1_list" --right "$R2_list" \
    --threads "$n_cores"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script transrate.sh"
date


#? TransRate paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4971766/
#? TransRate website: http://hibberdlab.com/transrate/metrics.html