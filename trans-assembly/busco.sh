#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --job-name=busco
#SBATCH --output=slurm-busco-%j.out


# PARSE COMMAND-LINE ARGS ------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run BUSCO to check a transcriptome or genome assembly."
  echo
  echo "Syntax: $0 -i <input-FASTA> -o <output-dir> -d <db-name> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Input FASTA file with the assembly"
  echo "    -o STRING         Output directory"
  echo "    -d STRING         Busco database name (see https://busco.ezlab.org/list_of_lineages.html)"
  echo
  echo "Other options:"
  echo "    -m STRING         Mode, i.e. assembly type               [default: 'transcriptome']"
  echo "                      Valid options: 'genome', 'transcripttome', or 'proteins'"
  echo "    -h                Print this help message and exit"
  echo
  echo "## Example:           $0 -i results/asssmbly/assembly.fa -o results/BUSCO -d bacteria_odb"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Busco docs: https://busco.ezlab.org/busco_userguide.html"
  echo
}

## Default parameter values
fa_in=""
outdir=""
busco_db=""
assembly_type=transcriptome

## Get parameter values
while getopts ':i:o:d:m:h' flag; do
    case "${flag}" in
    i) fa_in="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) busco_db="$OPTARG" ;;
    m) assembly_type="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Check input
[[ "$fa_in" = "" ]] && echo "## Please specify an input FASTA file with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## Please specify an output dir file with -o" >&2 && exit 1
[[ "$busco_db" = "" ]] && echo "## Please specify a Busco database name with -d" >&2 && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file $fa_in does not exist" >&2 && exit 1



# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/busco-env

## Bash strict mode
set -euo pipefail

## If needed, make input path absolute because we have to move into the outdir
[[ ! $fa_in =~ ^/ ]] && fa_in="$PWD"/"$fa_in"

## Get a sample/assembly ID from the filename
fileID=$(basename "$fa_in" | sed -E 's/.fn?as?t?a?//')

## Create output dir if needed
mkdir -p "$outdir"

## Report
echo
echo "## Starting script busco.sh"
date
echo
echo "## Input FASTA:               $fa_in"
echo "## Output dir:                $outdir"
echo "## BUSCO db:                  $busco_db"
echo "## Mode (assembly type):      $assembly_type"
echo "## Number of cores:           $SLURM_CPUS_PER_TASK"
echo "## Assembly ID (inferred):    $fileID"
echo -e "-------------------------------\n"


# RUN BUSCO -----------------------------------------------------------------
echo "## Now running busco..."

cd "$outdir" || exit 1

busco \
    -i "$fa_in" \
    -o "$fileID" \
    -l "$busco_db" \
    -m "$assembly_type" \
    -c "$SLURM_CPUS_PER_TASK" \
    --force


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh
echo -e "\n## Done with script busco.sh"
date

echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
