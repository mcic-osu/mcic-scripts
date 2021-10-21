#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --output=slurm-ASV_inference-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Process command-line arguments
config_file=$1

## Other variables/constants
n_cores=$SLURM_CPUS_PER_TASK

## Load modules
module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

## Run the R script
echo "## Config file: $config_file"
echo "## Number of cores: $n_cores"
echo -e "## Submitting script scripts/ASV_inference.R...\n"
Rscript mcic-scripts/metabarcoding/ASV_inference.R "$config_file" "$n_cores"
