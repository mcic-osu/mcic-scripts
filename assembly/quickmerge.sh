#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --output=slurm-quickmerge-%j.out
#SBATCH --mem=20G

## Load conda environment
module load python
conda activate /fs/project/PAS0471/jelmer/conda/quickmerge-env

## Command-line args
assembly1=$1
assembly2=$2
merged_assembly=$3

## Report
echo -e "\n## Starting script quickmerge.sh"
date
echo

## Create output dir
outdir=$(dirname "$merged_assembly")
mkdir -p "$outdir"

## Move into the outdir and change paths accordingly
assembly1="$PWD"/"$assembly1"
assembly2="$PWD"/"$assembly2"
merged_assembly=$(basename "$merged_assembly")
cd "$outdir" || exit

## Report
echo "## Assembly 1:           $assembly1"
echo "## Assembly 2:           $assembly2"
echo "## Merged assembly:      $merged_assembly"
echo

## Run Quickmerge
## https://github.com/mahulchak/quickmerge
merge_wrapper.py "$assembly1" "$assembly2" > "$merged_assembly"

echo -e "\n## Done with script"
date
