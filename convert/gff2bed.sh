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
  echo "## $0: Convert a GFF or GTF file to a BED file using 'bedops'."
  echo
  echo "## Syntax: $0 -i <input GFF/GTF> -o <output BED>..."
  echo
  echo "## Required options:"
  echo "## -i STRING        Input GFF/GTF file"
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
infile=""
bed=""

## Get command-line parameter values
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) infile="$OPTARG" ;;
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
[[ ! -f "$infile" ]] && echo "## ERROR: Input GFF/GTF file (-i) $infile does not exist" >&2 && exit 1
[[ $bed = "" ]] && echo "## ERROR: Please provide an output bed file with -o" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bedops-2.4.39

## Bash strict mode
set -euo pipefail

## Report
echo "## Input GFF file:                     $infile"
echo "## Output bed file:                    $bed"
echo

## Create output dir if needed
outdir=$(dirname "$bed")
mkdir -p "$outdir"


# CONVERT GFF TO bed -----------------------------------------------------------
if [[ $infile = *.gtf ]]; then
    echo "## Converting from GTF to BED..."
    gtf2bed < "$infile" > "$bed"
elif [[ $infile = *.gtf ]]; then
    echo "## Converting from GFF to BED"
    gff2bed < "$infile" > "$bed"
else
    echo "## ERROR: Not converting, make sure input file has extension '.gtf' or '.gff'"
fi


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing the input and output file:"
ls -lh "$gff" "$bed"
echo -e "\n## Done with script gff2bed.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
