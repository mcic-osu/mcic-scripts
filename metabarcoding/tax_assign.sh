#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --output=slurm-assign_tax-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Process command-line arguments
seqtab_rds=$1
taxa_rds=$2
algo=${3-decipher}

## Other variables
n_cores="$SLURM_CPUS_PER_TASK"

## Load modules
module load gnu/9.1.0
module load mkl/2019.0.5
module load R/4.0.2

## Checks
[[ ! -f $seqtab_rds ]] && echo "ERROR: Input file ($seqtab_rds) does not exist" && exit 1

## Report
echo "## Starting submission script dada2_qc_plots.sh..."
echo "## Seqtab RDS (input file): $seqtab_rds"
echo "## Taxa RDS (output file): $taxa_rds"
echo "## Taxonomic assignment algorithm: $algo"

## Run the R script
if [ "$algo" = "dada" ]; then
    echo -e "## Submitting script assign_tax_dada.R...\n"
    Rscript mcic-scripts/metabarcoding/assign_tax_dada.R "$seqtab_rds" "$taxa_rds" "$n_cores"
fi

if [ "$algo" = "decipher" ]; then
    echo -e "## Submitting script assign_tax_decipher.R...\n"
    Rscript mcic-scripts/metabarcoding/assign_tax_decipher.R "$seqtab_rds" "$taxa_rds" "$n_cores"
fi