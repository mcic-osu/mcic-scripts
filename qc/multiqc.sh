#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --out=slurm-multiqc-%j.out

# SETUP ------------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Run MultiQC."
  echo
  echo "## Syntax: $0 -i <input-dir> -o <output-dir> [-h]"
  echo
  echo "## Required options:"
  echo "## -i STRING     Input directory with e.g. FastQC output files"
  echo "## -o STRING     Output directory for MultiQC report"
  echo
  echo "## Other options:"
  echo "## -h            Print this help message and exit"
  echo
  echo "## Example command:"
  echo "## $0 -i results/fastqc -o results/multiqc"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
indir=""
outdir=""

## Parse command-line options
while getopts ':i:o:h' flag; do
  case "${flag}" in
  i) indir="$OPTARG" ;;
  o) outdir="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Report
echo "## Starting script multiqc.sh..."
date
echo

## Load software 
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/multiqc-1.12

## Bash strict mode
set -euo pipefail

## Input checks
[[ $indir = "" ]] && echo "## ERROR: Please specify input file with -i" >&2 && exit 1
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir does not exist" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify output dir with -i" >&2 && exit 1

## If necessary, create the output dir
mkdir -p "$outdir"

## Report
echo "## Input dir :       $indir"
echo "## Output dir:       $outdir"
echo -e "------------------\n"


# RUN MULTIQC ------------------------------------------------------------------
#? --interactive will ensure interactive plots, regardless of number of samples
#? --force will overwrite any old report

echo "## Starting MultiQC run..."
multiqc --interactive --force "$indir" -o "$outdir"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing output files:"
ls -lh "$outdir"
echo -e "\n## Done with script multiqc.sh"
date
echo
