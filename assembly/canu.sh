#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --job-name=canu
#SBATCH --output=slurm-canu-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run Canu to assemble a genome."
  echo
  echo "Syntax: $0 -i <FASTQ-file> -o <output-dir> -p <prefix> -s <genome-size> ..."
  echo
  echo "Required options:"
  echo "    -i STRING         Space-separated list of input FASTQ files"
  echo "    -p DIR            Custom prefix / name for the genome assembly"
  echo "    -o DIR            Output dir"
  echo "    -s STRING         Estimated genome size, e.g. '65m' for 65 Mbp"
  echo
  echo "Other options:"
  echo "    -a STRING         Time limit for Canu jobs, specify as HH:MM:SS     [default: 06:00:00]"
  echo "    -a STRING         Other argument(s) to pass to Canu"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i data/my.fastq.gz -o results/canu -p my_genome -s 250m"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Canu documentation: https://canu.readthedocs.io/en/latest/quick-start.html"
  echo
}

## Option defaults
infiles=""
prefix=""
outdir=""
genome_size=""
time_limit="6:00:00"
more_args=""

## Parse command-line options
while getopts ':i:o:p:t:a:s:h' flag; do
  case "${flag}" in
    i) infiles="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    t) time_limit="$OPTARG" ;;
    s) genome_size="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python
source activate /users/PAS0471/jelmer/miniconda3/envs/canu-env

## Bash script settings
set -euo pipefail

## Constants
CNS_MEMORY=4
GRID_MEMORY="--mem=20G"

## Arguments for Canu
slurm_options="--account=$SBATCH_ACCOUNT --time=$time_limit"
infile_arg="-nanopore ${infiles// / -nanopore }"

## Create the output directory
mkdir -p "$outdir"

## Report
echo "## Starting script canu.sh"
date
echo
echo "## Input FASTQ file(s):                  $infiles"
echo "## Prefix:                               $prefix"
echo "## Output dir:                           $outdir"
[[ $more_args != "" ]] && echo "## Other arguments to pass to Canu:    $more_args"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
echo "## Now running Canu..."
canu \
    -p "$prefix" \
    -d "$outdir" \
    genomeSize="$genome_size" \
    gridOptions="$slurm_options" \
    gridEngineMemoryOption="$GRID_MEMORY" \
    cnsMemory="$CNS_MEMORY" \
    $infile_arg $more_args


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script canu.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
