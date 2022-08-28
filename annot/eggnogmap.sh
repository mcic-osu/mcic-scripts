#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=60G
#SBATCH --job-name=eggnogmap
#SBATCH --output=slurm-eggnogmap-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run eggNOGmapper to functionally annotate a genome."
  echo
  echo "Syntax: $0 -i <protein-FASTA> -o <output-dir> -d <database-dir> ..."
  echo
  echo "Required options:"
  echo "    -i FILE           Input genome protein FASTA file"
  echo "    -d DIR            Pre-downloaded eggNOGmapper database dir"
  echo "    -o DIR            Output dir"
  echo
  echo "Other options:"
  echo "    -p STRING         Output prefix, such as an identifier for the genome   [default: 'genome']"
  echo "    -a STRING         Other argument(s) to pass to eggNOGmapper"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i proteins.fa -o results/eggnogmapper -d odb_prots.fa"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "eggNOGmapper documentation: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8"
  echo
}

## Option defaults
fa_in=""
db_dir=""
outdir=""
more_args=""
out_prefix=genome

## Parse command-line options
while getopts ':i:o:d:p:a:h' flag; do
  case "${flag}" in
    i) fa_in="$OPTARG" ;;
    d) db_dir="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    p) out_prefix="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Check input
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input file (-i) $fa_in does not exist" >&2 && exit 1
[[ ! -d "$db_dir" ]] && echo "## ERROR: Database dir (-d) $db_dir does not exist" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "## ERROR: Please specify an output dir with -o" >&2 && exit 1

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/eggnogg-env

## Bash strict mode
set -euo pipefail

## Other variables
temp_dir="$outdir"/tmp   # Alt: use compute node's $TMPDIR

## Make output dir
mkdir -p "$outdir" "$temp_dir"

## Report
echo "## Starting script eggnogmap.sh"
date
echo
echo "## Input FASTA file:                     $fa_in"
echo "## eggNOGmapper database dir:            $db_dir"
echo "## Output dir:                           $outdir"
echo "## Output prefix:                        $out_prefix"
[[ $more_args != "" ]] && echo "## Other arguments to pass to eggNOGmapper:    $more_args"
echo -e "--------------------\n"


# RUN EGGNOGMAPPER -------------------------------------------------------------
echo "## Now running eggNOGmapper..."
emapper.py \
    -i "$fa_in" \
    --data_dir "$db_dir" \
    --output_dir "$outdir" \
    --output "$out_prefix" \
    -m diamond \
    --go_evidence all \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --temp_dir "$temp_dir" \
    --override $more_args

#? Non-default options
# --go_evidence all => default is to use only non-electronic terms (`non-electronic`), see https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8

## Other options
# --pfam_realign denovo \ #! Needs some HMMer server setup


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo
echo "## Done with script eggnogmap.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
