#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-assign_tax-%j.out

## Process command-line arguments
seqtab_rds=$1
taxa_rds=$2
algo=${3-dada}

## Other variables
SCRIPT_DADA=mcic-scripts/metabarcoding/tax_assign_dada.R
SCRIPT_DECI=mcic-scripts/metabarcoding/tax_assign_deci.R
n_cores="$SLURM_CPUS_PER_TASK"

## Load modules
module load R/4.0.2-gnu9.1

## Checks
[[ ! -f $seqtab_rds ]] && echo "ERROR: Input file ($seqtab_rds) does not exist" && exit 1

## Report
echo "## Starting submission script dada2_qc_plots.sh..."
echo "## Seqtab RDS (input file): $seqtab_rds"
echo "## Taxa RDS (output file): $taxa_rds"
echo "## Taxonomic assignment algorithm: $algo"

## Run the R script
if [ "$algo" = "dada" ]; then
    echo -e "## Running script assign_tax_dada.R...\n"
    Rscript "$SCRIPT_DADA" "$seqtab_rds" "$taxa_rds" "$n_cores"
fi

if [ "$algo" = "decipher" ]; then
    echo -e "## Running script assign_tax_decipher.R...\n"
    Rscript "$SCRIPT_DECI" "$seqtab_rds" "$taxa_rds" "$n_cores"
fi

echo "## Done with script tax_assign.sh"
date
echo
