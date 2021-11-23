#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --output=slurm-fastqc-%j.out

## Bash strict settings
set -ueo pipefail  # Bash strict settings

## Load software
module load fastqc

## Help function
Help() {
  echo
  echo "## $0: Run FastQC for a single FASTQ file."
  echo
  echo "## Syntax: $0 -i <input-FASTQ> -o <output-dir> [-h]"
  echo "## Options:"
  echo "## -h       Print this help message"
  echo "## -i STR   Input FASTQ file (REQUIRED)"
  echo "## -o STR   Output directory (REQUIRED)"
  echo "## Example: $0 -i data/fastq/my.fastq.gz -o results/fastqc"
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
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Input checks
[[ $infile = "" ]] && echo "## ERROR: Please specify input file with -i" && exit 1
[[ ! -f "$infile" ]] && echo "## ERROR: Input file does not exist" && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify output dir with -i" && exit 1

## Create the output directory if it doesn't already exist
mkdir -p "$outdir"

## Report
echo "## Starting script fastqc.sh..."
date
echo "Input FASTQ file:       $infile"
echo "Output dir:             $outdir"
echo -e "--------------------\n"

## Run FastQC
fastqc --outdir="$outdir" "$infile"

## Report
echo -e "\n## Listing output files:"
ls -lh "$outdir"
echo -e "\n## Done with script fastqc.sh"
date