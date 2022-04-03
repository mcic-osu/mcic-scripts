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
  echo "## $0: Run BUSCO to check a transcriptome or genome assembly."
  echo
  echo "## Syntax: $0 -i <input-FASTA> -o <output-dir> -d <db-name> ..."
  echo
  echo "## Required options:"
  echo "## -i STRING     Input FASTA file with the assembly"
  echo "## -o STRING     Output directory"
  echo "## -d STRING     Busco database name (see https://busco.ezlab.org/list_of_lineages.html)"
  echo
  echo "## Other options:"
  echo "## -t STRING     Assembly type                  [default: 'transcriptome']"
  echo "## -h            Print this help message and exit"
  echo
  echo "## Example: $0 -i results/asssmbly/assembly.fa -o results/BUSCO -d bacteria_odb"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameter values
fa_in=""
outdir=""
busco_db=""
type=transcriptome

## Get parameter values
while getopts ':i:o:d:t:h' flag; do
    case "${flag}" in
    i) fa_in="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    d) busco_db="$OPTARG" ;;
    t) type="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo "## Starting script busco.sh"
date
echo

## Check input
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file $fa_in does not exist" && exit 1


# SOFTWARE ---------------------------------------------------------------------
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/busco-env


# OTHER SETUP ------------------------------------------------------------------
## Bash strict mode
set -euo pipefail

## Other paramaters
n_cores=$SLURM_CPUS_PER_TASK

## Report
echo "## Input FASTA:           $fa_in"
echo "## Output dir:            $outdir"
echo "## BUSCO db:              $busco_db"
echo "## Assembly type:         $type"
echo "## Number of cores:       $n_cores"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN BUSCO -----------------------------------------------------------------
echo "## Now running busco..."
busco \
    -i "$fa_in" \
    -o "$outdir" \
    -l "$busco_db" \
    -m "$type" \
    -c "$n_cores" \
    --force


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script busco.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo


# INFO -------------------------------------------------------------------------
#? User guide: https://busco.ezlab.org/busco_userguide.html
#? Conda: https://anaconda.org/bioconda/busco