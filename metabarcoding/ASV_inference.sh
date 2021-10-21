#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --output=slurm-ASV_inference-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Process command-line arguments
fastq_indir=$1
outdir=$2
config_file=$3

## Other variables/constants
n_cores=$SLURM_CPUS_PER_TASK

## Load modules
module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

## Checks
[[ ! -d $fastq_indir ]] && echo "ERROR: FASTQ input dir ($fastq_indir) does not exist" && exit 1
[[ ! -f $config_file ]] && echo "ERROR: Config file ($config_file) does not exist" && exit 1

## Report
echo "## Starting submission script ASV_inference.sh..."
echo "## FASTQ input dir: $fastq_indir"
echo "## Output dir: $outdir"
echo "## Config file: $config_file"
echo "## Number of cores: $n_cores"

## Run the R script
echo -e "## Submitting script scripts/ASV_inference.R...\n"
Rscript mcic-scripts/metabarcoding/ASV_inference.R \
    "$fastq_indir" "$outdir" "$config_file" "$n_cores"
