#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --job-name=filtlong
#SBATCH --output=slurm-filtlong-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run FiltLong to filter long reads."
  echo
  echo "Syntax: $0 -i <FASTQ-file> -o <output-dir> -t <nr-bp> ..."
  echo
  echo "Required options:"
  echo "    -i FILE           Input FASTQ file"
  echo "    -o DIR            Output dir"
  echo "    -t STRING         Target number of bases"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to FiltLong"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i data/my.fastq.gz -o results/filtlong -t 1000000"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "FiltLong documentation: https://github.com/rrwick/Filtlong"
  echo
}

## Option defaults
infile=""
outdir=""
target_bases=""
more_args=""

## Parse command-line options
while getopts ':i:o:t:a:h' flag; do
  case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    t) target_bases="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python
conda activate /fs/project/PAS0471/jelmer/conda/filtlong-env

## Bach script settings
set -euo pipefail

## Check if input dir isn't the same as the output dir
indir=$(dirname "$infile")
[[ "$indir" = "$outdir" ]] && echo "## ERROR: Input dir can't be the same as output dir" && exit 1

## Determine output file
outfile="$outdir"/$(basename "$infile")

## Create the output directory
mkdir -p "$outdir"

## Report
echo "## Starting script filtlong.sh"
date
echo
echo "## Input FASTQ file:                     $infile"
echo "## Output dir:                           $outdir"
echo "## Target nr of bases:                   $target_bases"
[[ $more_args != "" ]] && echo "## Other arguments to pass to FiltLong:    $more_args"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
filtlong \
    --target_bases "$target_bases" \
    "$infile" |
    gzip > "$outfile"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing the output file:"
ls -lh "$outfile"
echo -e "\n## Done with script filtlong.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
