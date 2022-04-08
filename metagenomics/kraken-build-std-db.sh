#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=20
#SBATCH --output=slurm-kraken-build-std-db-%j.out

## Help function
Help() {
    echo
    echo "## $0: Download and build the standard Kraken db."
    echo
    echo "## Syntax: $0 -d <kraken-db-dir> ..."
    echo 
    echo "## Required options:"
    echo "## -d     Directory for the Kraken db"
    echo
    echo "## Other options"
    echo "## -h     Print this help message and exit"
    echo
    echo "## Example: $0 -d refdata/kraken/my-db"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}


# SET-UP -----------------------------------------------------------------------
## Report
echo -e "\n## Starting script kraken-build-std-db.sh..."
date

## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/kraken2-env

## Bash strict mode
set -euo pipefail

## Option defaults
db_dir=""

## Get command-line options
while getopts ':d:h' flag; do
    case "${flag}" in
    d) db_dir="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Check input
[[ "$db_dir" = "" ]] && echo -e "\n## ERROR: must specify a db dir with -d" && exit 1

## Make DB directory
mkdir -p "$db_dir"

## Report
echo
echo "## Database dir:                 $db_dir"
echo -e "------------------\n"


# BUILD THE KRAKEN2 DATABASE ---------------------------------------------------
echo "## Building the Kraken database..."
kraken2-build --standard \
    --db "$db_dir" \
    --threads "$SLURM_CPUS_ON_NODE"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing file in DB dir:"
ls -lh "$db_dir"
echo -e "\n## Done with script kraken-build-std-db.sh"
date