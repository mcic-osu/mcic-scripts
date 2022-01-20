#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=kallisto
#SBATCH --output=slurm-kallisto-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Quantify transcripts with Kallisto."
  echo
  echo "## Syntax: $0 -i <R1-FASTQ-infile> -o <BAM-outdir> -r <ref-index-dir> [ -a <gff-file> ... ] [-h]"
  echo
  echo "## Options:"
  echo "##    -h        Print this help message"
  echo "##    -i STR    R1 FASTQ input file (REQUIRED; note that the name of the R2 file will be inferred by the script.)"
  echo "##    -r STR    Reference index (REQUIRED)"
  echo "##    -o STR    Output dir (REQUIRED)"
  echo
  echo "## Example: $0 -i data/fastq/S01_L001_R1.fastq.gz -o results/kallisto -r refdata/trans.idx"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Parse command-line options
while getopts ':i:o:r::h' flag; do
  case "${flag}" in
  i) R1="$OPTARG" ;;
  o) outdir="$OPTARG" ;;
  r) index="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done


# SETUP ---------------------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/kallisto-env

## Bash strict mode
set -euo pipefail

## Process parameters
R2=${R1/_R1_/_R2_}                               # R2 FASTQ file
sample_id=$(basename "$R1" | sed 's/_R1_.*//')   # Sample ID - for output files
outdir_full=$outdir/"$sample_id"

## Report
echo "## Starting script kallisto_quant.sh"
date
echo
echo "## Input R1 FASTQ file:                          $R1"
echo "## Transcriptome index:                          $index"
echo "## Output dir - base:                            $outdir"
echo
echo "## Output dir - full:                            $outdir_full"
echo "## Sample ID (as inferred by the script):        $sample_id"
echo "## R2 FASTQ file (as inferred by the script):    $R2"
echo -e "------------------------\n"

## Check inputs
[[ ! -f "$R1" ]] && echo "## ERROR: Input file R1_in (-i) $R1 does not exist" >&2 && exit 1
[[ ! -f "$R2" ]] && echo "## ERROR: Input file R2_in $R2 does not exist" >&2 && exit 1

## Create output dir if needed
mkdir -p "$outdir"


# QUANTIFY ---------------------------------------------------------------------
kallisto quant \
    -i "$index" \
    -o "$outdir_full" \
    -b 100 \
    "$R1" "$R2"

#? -b = number of bootstraps

# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"

echo -e "\n## Listing output files:"
ls -lh "$outdir_full"

echo -e "\n## Done with script kallisto_quant.sh"
date
