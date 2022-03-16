#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --output=slurm-fastqc-%j.out


# SETUP ------------------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Run FastQC for a single FASTQ file."
  echo
  echo "## Syntax: $0 -i <input-FASTQ> -o <output-dir> [-h]"
  echo
  echo "## Required options:"
  echo "## -i STR     Input FASTQ file"
  echo "## -o STR     Output directory"
  echo
  echo "## Other options:"
  echo "## -h         Print this help message"
  echo
  echo "## Example command:"
  echo "## $0 -i data/fastq/my.fastq.gz -o results/fastqc"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Option defaults
infile=""
outdir=""

## Parse command-line options
while getopts ':i:o:h' flag; do
  case "${flag}" in
  i) infile="$OPTARG" ;;
  o) outdir="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Report
echo "## Starting script fastqc.sh..."
date
echo

## Load software
module load fastqc

## Bash strict settings
set -euo pipefail

## Input checks
[[ $infile = "" ]] && echo "## ERROR: Please specify input file with -i" >&2 && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input file does not exist" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify output dir with -o"  >&2 && exit 1

## Create the output directory if it doesn't already exist
mkdir -p "$outdir"

## Report
echo "## Input FASTQ file:       $infile"
echo "## Output dir:             $outdir"
echo -e "--------------------\n"


# RUN FASTQC -------------------------------------------------------------------
fastqc --outdir="$outdir" "$infile"


# WRAP UP ----------------------------------------------------------------------
sample_id=$(basename "$infile" | sed 's/.fastq.*//')
echo -e "\n## Listing output files:"
ls -lh "$outdir"/"$sample_id"*fastqc*

echo -e "\n## Done with script fastqc.sh"
date
echo
