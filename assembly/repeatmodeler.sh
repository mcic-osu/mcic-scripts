#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --job-name=repeatmodeler
#SBATCH --output=slurm-repeatmodeler-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run RepeatModeler to search for repeats in a genome."
  echo
  echo "Syntax: $0 -i <genome-FASTA> -o <output-dir>..."
  echo
  echo "Required options:"
  echo "    -i FILE           Genome (nucleotide) FASTA file"
  echo "    -o DIR            Output dir"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to RepeatModeler"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i my_genome.fa -o results/repeatmodeler"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "RepeatModeler documentation: https://github.com/Dfam-consortium/RepeatModeler"
  echo "RepeatModeler paper: https://www.pnas.org/doi/10.1073/pnas.1921046117"
  echo
}

## Option defaults
genome_fa=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
  case "${flag}" in
    i) genome_fa="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

# SETUP ------------------------------------------------------------------------
## Check input
[[ ! -f "$genome_fa" ]] && echo "## ERROR: Input file (-i) $genome_fa does not exist" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please provide an output dir with -o" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/repeatmodeler-2.0.3

## Bash script settings
set -euo pipefail

## Get genome ID for database prefix
genomeID=$(basename "$genome_fa" | sed -E 's/\.fn?as?t?a?//')

## Report
echo
echo "## Starting script repeatmodeler.sh"
date
echo
echo "## Input file (genome FASTA):                  $genome_fa"
echo "## Output dir:                                 $outdir"
[[ $more_args != "" ]] && echo "## Other arguments to pass to RepeatModeler:    $more_args"
echo -e "--------------------\n"

## Make output dir
mkdir -p "$outdir"

# RUN --------------------------------------------------------------------------
echo "## Now building the RepeatModeler database..."
BuildDatabase \
    -name "$outdir"/"$genomeID" \
    "$genome_fa"

echo -e "\n## Now runnning RepeatModeler..."
RepeatModeler \
    -database "$outdir"/"$genomeID" \
    $more_args "$genome_fa"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script repeatmodeler.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
