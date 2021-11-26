#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --output=slurm-make_ps-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Process command-line arguments
seqtab_rds=$1      # Input seqtab object RDS file
taxa_rds=$2        # Input taxa object RDS file
tree_rds=$3        # Input tree object RDS file
sampledata_file=$4 # Sample data text file
ps_rds=$5          # Output phyloseq object RDS file

## Other variables
PS_SCRIPT=mcic-scripts/metabarcoding/ps_make.R
n_cores="$SLURM_CPUS_PER_TASK"

## Load modules
module load R/4.0.2-gnu9.1

## Checks
[[ ! -f $seqtab_rds ]] && echo "ERROR: Input file ($seqtab_rds) does not exist" && exit 1
[[ ! -f $taxa_rds ]] && echo "ERROR: Input file ($taxa_rds) does not exist" && exit 1
[[ ! -f $tree_rds ]] && echo "ERROR: Input file ($tree_rds) does not exist" && exit 1
[[ ! -f $sampledata_file ]] && echo "ERROR: Input file ($sampledata_file) does not exist" && exit 1

## Report
echo "## Starting submission script ps_make.sh..."
echo "## Seqtab RDS file (input): $seqtab_rds"
echo "## Taxa RDS file (input): $taxa_rds"
echo "## Tree RDS file (input): $tree_rds"
echo "## Sample metadat file (input): $sampledata_file"
echo "## Phyloseq RDS file (output): $ps_rds"
echo

## Run the R script
Rscript "$PS_SCRIPT" "$seqtab_rds" "$taxa_rds" "$tree_rds" "$sampledata_file" "$ps_rds" "$n_cores"
