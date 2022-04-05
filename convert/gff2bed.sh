#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --job-name=gff2bed
#SBATCH --output=slurm-gff2bed-%j.out


# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Convert a GFF file to a BED file using 'bedops'."
  echo
  echo "## Syntax: $0 -i <input GFF> -o <output BED>..."
  echo
  echo "## Required options:"
  echo "## -i STRING        Input GFF file"
  echo "## -o STRING        Output BED file"
  echo
  echo "## Other options:"
  echo "## -h               Print this help message and exit"
  echo
  echo "## Example: $0 -i my.gff -o my.bed"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameter values
gff=""
bed=""

## Get command-line parameter values
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) gff="$OPTARG" ;;
    o) bed="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Report
echo -e "\n## Starting script gff2bed.sh"
date
echo

## Test parameter values
[[ ! -f "$gff" ]] && echo "## ERROR: Input GFF file (-i) $gff does not exist" >&2 && exit 1
[[ $bed = "" ]] && echo "## ERROR: Please provide an output bed file with -o" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bedops-2.4.39

## Bash strict mode
set -euo pipefail

## Report
echo "## Input GFF file:                     $gff"
echo "## Output bed file:                    $bed"
echo

## Create output dir if needed
outdir=$(dirname "$bed")
mkdir -p "$outdir"


# CONVERT GFF TO bed -----------------------------------------------------------
gff2bed < "$gff" > "$bed"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing the input and output file:"
ls -lh "$gff" "$bed"
echo -e "\n## Done with script gff2bed.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
