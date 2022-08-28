#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=20
#SBATCH --job-name=orthofinder
#SBATCH --output=slurm-orthofinder-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run OrthoFinder to find orthologs between genomes."
  echo
  echo "Syntax: $0 -i <input-dir> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i DIR            Input dir with protein FASTA files, one per genome"
  echo "    -o DIR            Output dir"
  echo "                      WARNING: If output dir exists, it will be emptied prior to running OrthoFinder"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to OrthoFinder"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i data/genomes -o results/orthofinder"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "OrthoFinder documentation: https://github.com/davidemms/OrthoFinder and https://davidemms.github.io/"
  echo "OrthoFinder paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y"
  echo
}

## Option defaults
indir=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
  case "${flag}" in
    i) indir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Check input
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir (-i) $indir does not exist" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/orthofinder-2.5.4

## Bash strict mode
set -euo pipefail

## Remove output dir
rm -r "$outdir"

## Report
echo
echo "## Starting script orthofinder.sh"
date
echo
echo "## Input dir:                    $indir"
echo "## Output dir:                   $outdir"
[[ $more_args != "" ]] && echo "## Other arguments for OrthoFinder:    $more_args"
echo -e "--------------------\n"

# RUN --------------------------------------------------------------------------
echo "## Starting Orthofinder run..."
orthofinder \
    -f "$indir" \
    -o "$outdir" \
    -t "$SLURM_CPUS_PER_TASK" \
    -a "$SLURM_CPUS_PER_TASK" $more_args


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script orthofinder.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
