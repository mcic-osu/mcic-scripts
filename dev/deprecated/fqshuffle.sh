#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --time=30
#SBATCH --account=PAS0471
#SBATCH --job-name=fqshuffle
#SBATCH --output=slurm-fqshuffle-%j.out


# SET-UP & PARSE ARGS ----------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/bbmap-38.96

## Bash strict settings
set -euo pipefail

## Help function
Help() {
    echo "## fqshuffle.sh: script to shuffle FASTQ files using BBmap shuffle.sh"
    echo
    echo "## Syntax: fqshuffle.sh -i <R1_in> -I <R2_in> -o <R1_out> -o <R2_out> [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -i     R1 input file (REQUIRED)"
    echo "## -I     R2 input file (REQUIRED)"
    echo "## -o     R1 output file (REQUIRED)"
    echo "## -O     R2 output file (REQUIRED)"
    echo
}

## Parse command-line options
while getopts ':i:I:o:O:h' flag; do
  case "${flag}" in
    i)  R1_in="$OPTARG" ;;
    I)  R2_in="$OPTARG" ;;
    o)  R1_out="$OPTARG" ;;
    O)  R2_out="$OPTARG" ;;
    h)  Help && exit 0 ;;
	\?) echo "## ERROR: Invalid option" >&2 && exit 1 ;;
	:)  echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Check input - error out if neither n_reads or prop_reads is provided
[[ ! -f "$R1_in" ]] && echo "ERROR: Input file $R1_in does not exist" >&2 && exit 1
[[ ! -f "$R2_in" ]] && echo "ERROR: Input file $R2_in does not exist" >&2 && exit 1

## Create output dir if needed
outdir=$(dirname "$R1_in")
mkdir -p "$outdir"

## Report
echo -e "\n## Starting script fqshuffle.sh..."
date
echo "## Input R1:       $R1_in"
echo "## Input R2:       $R2_in"
echo "## Output R1:      $R1_out"
echo "## Output R2:      $R2_out"
echo
echo "## Listing input files:"
ls -lh "$R1_in"
ls -lh "$R2_in"
echo -e "-------------------------\n"


# RUN SEQTK TO SUBSAMPLE FASTQ -------------------------------------------------
shuffle.sh \
    in="$R1_in" in2="$R2_in" \
    out="$R1_out" out2="$R2_out"


# REPORT AND FINALIZE --------------------------------------------------------
echo -e "\n----------------------------------"
echo -e "## Listing output files:"
ls -lh "$R1_out"
ls -lh "$R2_out"

echo -e "\n## Done with script."
date
