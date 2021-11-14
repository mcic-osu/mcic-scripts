#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-dada2-qc-plots-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

## Process command-line arguments
nreads_file=$1
outdir=$2

## Other scripts
QC_SCRIPT=mcic-scripts/metabarcoding/ASV_qcplots.R

## Load modules
module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

## Checks
[[ ! -f $nreads_file ]] && echo "ERROR: Input file ($nreads_file) does not exist" && exit 1

## Report
echo "## Starting submission script ASV_qcplots.sh..."
echo "## Input file: $nreads_file"
echo "## Output dir: $outdir"

## Run the R script
Rscript "$QC_SCRIPT" "$nreads_file" "$outdir"
