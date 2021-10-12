#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --output=slurm-ASV_inference-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Process command-line arguments
fastq_indir=$1    # Dir with input FASTQ files
outdir=$2         # Dir for output 
n_samples=$3      # Number of samples to run ("all" => all samples, otherwise specify an integer)
trunc_f=${4-180}  # Truncate forward reads at this length
trunc_r=${4-180}  # Truncate reverse reads at this length

## Other variables/constants
n_cores=$SLURM_CPUS_PER_TASK

## Load modules
module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

## Run the R script
Rscript scripts/ASV_inference.R \
    "$fastq_indir" "$outdir" "$n_cores" "$n_samples" "$trunc_f" "$trunc_r"
