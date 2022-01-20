#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=kallisto_index
#SBATCH --output=slurm-kallisto_index-%j.out


# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "## $0: Index a transcriptome with Kallisto."
  echo
  echo "## Syntax: $0 -i <transcriptome-fasta> -o <transcriptome-index> [-h]"
  echo
  echo "## Options:"
  echo "##    -h        Print this help message"
  echo "##    -i STR    Input: transcriptomoe FASTA file (REQUIRED)"
  echo "##    -o STR    Output: transcriptome index file (REQUIRED)"
  echo
  echo "## Example: $0 -i -o results/trinity/trans.fa -r results/kallisto/trans.idx"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
}

## Parse command-line options
while getopts ':i:o:r::h' flag; do
  case "${flag}" in
  i) transcriptome="$OPTARG" ;;
  o) index="$OPTARG" ;;
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

## Report
echo "## Starting script kallisto_index.sh"
date
echo
echo "## Transcriptome FASTA (input):                  $transcriptome"
echo "## Transcriptome index (output):                 $index"
echo -e "------------------------\n"

## Check inputs
[[ ! -f "$transcriptome" ]] && echo "## ERROR: Input transcriptome file (-i) $transcriptome does not exist" >&2 && exit 1

## Create output dirs if needed
outdir=$(dirname "$index")
mkdir -p "$outdir"


# INDEX ------------------------------------------------------------------------
kallisto index -i "$index" "$transcriptome"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"

echo -e "\n## Listing output file:"
ls -lh "$index"

echo -e "\n## Done with script kallisto_index.sh"
date
