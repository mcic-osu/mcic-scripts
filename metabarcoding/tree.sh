#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --output=slurm-tree-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Process command-line arguments
seqtab_rds=$1
tree_rds=$2

## Other variables/constants
n_cores=$SLURM_CPUS_PER_TASK

## Load modules
module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

## Checks
[[ ! -f $seqtab_rds ]] && echo "ERROR: Input file ($seqtab_rds) does not exist" && exit 1

## Report
echo "## Starting submission script tree.sh..."
echo "## Sequence table RDS file (input): $seqtab_rds"
echo "## Tree RDS file (output): $tree_rds"
echo "## Number of cores: $n_cores"

## Run the R script
echo -e "## Submitting script tree.R...\n"
Rscript mcic-scripts/metabarcoding/tree.R "$seqtab_rds" "$tree_rds" "$n_cores"