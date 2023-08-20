#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --job-name=gff2gtf
#SBATCH --output=slurm-gff2gtf-%j.out


# PARSE ARGS -------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Convert a GFF file to a GTF file using 'gffread'."
  echo
  echo "## Syntax: $0 -i <input GFF> -o <output GTF>..."
  echo
  echo "## Required options:"
  echo "## -i STRING        Input GFF file"
  echo "## -o STRING        Output GTF file"
  echo
  echo "## Other options:"
  echo "## -h               Print this help message and exit"
  echo
  echo "## Example: $0 -i my.gff -o my.gtf"
  echo "## To submit to the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Default parameter values
gff=""
gtf=""

## Get command-line parameter values
while getopts ':i:o:h' flag; do
    case "${flag}" in
    i) gff="$OPTARG" ;;
    o) gtf="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done


# OTHER SETUP ------------------------------------------------------------------
## Report
echo -e "\n## Starting script gff2gtf.sh"
date
echo

## Test parameter values
[[ ! -f "$gff" ]] && echo "## ERROR: Input GFF file (-i) $gff does not exist" >&2 && exit 1
[[ $gtf = "" ]] && echo "## ERROR: Please provide an output GTF file with -o" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/gffread-0.12.7

## Bash strict mode
set -euo pipefail

## Report
echo "## Input GFF file:                     $gff"
echo "## Output GTF file:                    $gtf"
echo

## Create output dir if needed
outdir=$(dirname "$gtf")
mkdir -p "$outdir"


# CONVERT GFF TO GTF -----------------------------------------------------------
gffread "$gff" -T -o "$gtf"

#? Command from nf-core RNAseq:
#gffread GCA_003693625.1.gff --keep-exon-attrs -F -T -o GCA_003693625.1.gtf

# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing the input and output file:"
ls -lh "$gff" "$gtf"
echo -e "\n## Done with script gff2gtf.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
