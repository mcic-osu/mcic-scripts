#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-assign_tax-%j.out

## Help function
Help() {
    echo
    echo "## $0: Assign taxonomy to a dada-derived sequence table."
    echo
    echo "## Syntax: $0 -i <input-R1-FASTQ> -o <output-dir> ..."
    echo
    echo "## Required options:"
    echo "## -i STR     Input sequence table RDS file"
    echo "## -o STR     Output taxonomic table RDS file"
    echo
    echo "## Other options:"
    echo "## -a STR     Algorithm to use: either 'dada' or 'deci' (default: 'dada'; 'deci' is DECIPHER's IdTaxa function)"
    echo "## -h         Print this help message"
    echo
    echo "## Example: $0 -i results/dada/seqtab.rds -o results/taxonomy/tax.rds -a deci"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}


# SETUP ------------------------------------------------------------------------
## Report
echo "## Starting script tax_assign.sh..."
date
echo

## Load software
module load R/4.0.2-gnu9.1

## Bash strict settings
set -euo pipefail

## Option defaults
seqtab_rds=""
taxa_rds=""
algo="dada"

## Parse command-line options
while getopts ':i:o:a:h' flag; do
    case "${flag}" in
    i) seqtab_rds="$OPTARG" ;;
    o) taxa_rds="$OPTARG" ;;
    a) algo="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Other variables
SCRIPT_DADA=mcic-scripts/metabarcoding/tax_assign_dada.R
SCRIPT_DECI=mcic-scripts/metabarcoding/tax_assign_deci.R
n_cores="$SLURM_CPUS_PER_TASK"

## Check input
[[ ! -f $seqtab_rds ]] && echo "## ERROR: Input file ($seqtab_rds) does not exist" >&2 && exit 1
[[ $algo != "dada" && $algo != "decipher" ]] && echo "## ERROR: algorith should be 'dada' or 'decipher', not $algo" >&2 && exit 1

## Report
echo "## Seqtab RDS (input file):          $seqtab_rds"
echo "## Taxa RDS (output file):           $taxa_rds"
echo "## Taxonomic assignment algorithm:   $algo"
echo


# RUN ONE OF THE ASSIGNMENT SCRIPTS --------------------------------------------
## Run the R script
if [ "$algo" = "dada" ]; then
    echo -e "## Running script assign_tax_dada.R...\n"
    Rscript "$SCRIPT_DADA" -i "$seqtab_rds" -o "$taxa_rds" -c "$n_cores"
fi

if [ "$algo" = "decipher" ]; then
    echo -e "## Running script assign_tax_decipher.R...\n"
    Rscript "$SCRIPT_DECI" -i "$seqtab_rds" -o "$taxa_rds" -c "$n_cores"
fi

## Report
echo -e "\n## Done with script tax_assign.sh"
date
echo
