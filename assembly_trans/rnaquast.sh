#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --job-name=rnaquast
#SBATCH --output=slurm-rnaquast-%j.out


# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run rnaQUAST to evaluate transcriptome assembly quality."
  echo
  echo "Syntax: $0 -i <input-FASTA> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "-i STRING         One or more input FASTA files (transcriptome assemblies)"
  echo "                     When providing more than 1 assembly, use a space-separated list"
  echo "-o STRING         Output dir"
  echo
  echo "Other options:"
  echo "-x STRING         Extension of assembly files          [default: fasta]"
  echo "-a STRING         Other argument(s) to pass to CheckM"
  echo "-h                Print this help message and exit"
  echo
  echo "Example: $0 -i results/assembly -o results/checkm -x fa"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "rnaQUAST documentation: https://github.com/ablab/rnaquast"
  echo
}

## Option defaults
infiles=""
outdir=""
more_args=""

## Parse command-line options
while getopts ':i:o:a:h' flag; do
  case "${flag}" in
    i) infiles="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done

## Other parameters
n_cores=$SLURM_CPUS_PER_TASK

## Check input
[[ "$infiles" = "" ]] && echo "## ERROR: Please specify one or more input files with -i" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/rnaquast-env

## Bash strict mode
set -euo pipefail

## Report
echo "## Starting script rnaquast.sh"
date
echo "## Input FASTA file(s):    $infiles"
echo "## Output dir:             $outdir"
echo "## Number of cores:        $n_cores"
echo -e "------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# RUN RNAQUAST -----------------------------------------------------------------
echo "## Now running rnaQUAST..."
rnaQUAST.py \
    --transcripts $infiles \
    --output_dir "$outdir" \
    --threads "$n_cores" \
    --strand_specific \
    --gene_mark \
    $more_args

#? Not running BUSCO with rnaQUAST because it failed


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script rnaquast.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
