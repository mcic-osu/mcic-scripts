#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --job-name=interproscan
#SBATCH --output=slurm-interproscan-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run InterProScan to functionally annotate a genome."
  echo
  echo "Syntax: $0 -i <protein-FASTA> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Protein FASTA file"
  echo "    -o STRING         Output dir"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to InterProScan"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i my_genome.fa -o results/braker -d odb_prots.fa"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "InterProScan documentation: https://interproscan-docs.readthedocs.io/"
  echo
}

## Option defaults
protein_fa=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
  case "${flag}" in
    i) protein_fa="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Check input
[[ ! -f "$protein_fa" ]] && echo "## ERROR: Protein file (-d) $protein_fa does not exist" >&2 && exit 1

## Braker2 conda env which contains everything except GeneMark-EX and ProtHint 
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/interproscan-5.55

## Bash strict mode
set -euo pipefail

## Make output dir
mkdir -p "$outdir"

## Report
echo "## Starting script interproscan.sh"
date
echo
echo "## Protein FASTA file:                   $protein_fa"
echo "## Output dir:                           $outdir"
[[ $more_args != "" ]] && echo "## Other arguments to pass to InterProScan:    $more_args"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
echo "## Now running InterProScan..."
interproscan.sh \
    -i "$protein_fa" \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --output-dir "$outdir" \
    --goterms \
    --tempdir "$TMPDIR" \
    --disable-precalc

# --pathways -- add this option?

# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script interproscan.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo


# DEV NOTES --------------------------------------------------------------------
## Setup
# python3 /fs/project/PAS0471/jelmer/conda/interproscan-5.55/share/InterProScan/initial_setup.py
# Note: above didn't work, but running interproscan.sh without args seemed to do the setup!
