#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-entap-config-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Configure EnTAP before transcriptome assembly."
  echo
  echo "Syntax: $0 -i <genome-FASTA> -c <config-file> -d <db-file> -o <output-dir>..."
  echo
  echo "Required options:"
  echo "-c STRING         Config file"
  echo "-d STRING         Space-separated list of reference protein FASTA files"
  echo "-o STRING         Output dir"
  echo
  echo "Other options:"
  echo "-a STRING         Other argument(s) to pass to EnTAP"
  echo "-h                Print this help message and exit"
  echo
  echo "Example:    $0 -c entap_config.ini -d uniprot_sprot.fa -o results/entap"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "EnTAP documentation: https://entap.readthedocs.io/en/latest/"
  echo "EnTAP paper: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13106"
  echo
}

## Option defaults
dbs=""
outdir=""
config_file=""
more_args=""

## Parse command-line options
while getopts ':d:o:c:a:h' flag; do
  case "${flag}" in
    c) config_file="$OPTARG" ;;
    d) dbs="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/entap-0.10.8

## Bash strict mode
set -euo pipefail

## Build database arg
db_arg="-d ${dbs// / -d }"

## Report
echo
echo "## Starting script entap_config.sh"
date
echo
echo "## Input database FASTAs:                $dbs"
echo "## Input database FASTA arg:             $db_arg"
echo "## Output database dir:                  $outdir"
echo "## Config file:                          $config_file"
echo "## Other arguments to pass to EnTAP:     $more_args"
echo -e "--------------------\n"


# RUN ENTAP CONFIG -------------------------------------------------------------
echo "## Now running EnTAP config...."
EnTAP --config \
    --ini "$config_file" \
    --out-dir "$outdir" \
    -t "$SLURM_CPUS_PER_TASK" \
    $db_arg $more_args


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script entap_config.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
