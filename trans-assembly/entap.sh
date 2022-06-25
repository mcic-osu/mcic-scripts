#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-entap-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run EnTAP to annotate a transcriptome assembly."
  echo
  echo "Syntax: $0 -i <genome-FASTA> -c <config-file> -d <db-file> -o <output-dir>..."
  echo
  echo "Required options:"
  echo "-i STRING         Input assembly (nucleotide) FASTA file"
  echo "-c STRING         Config file"
  echo "-o STRING         Output dir"
  echo "-d STRING         Space-separated list of DIAMOND database files ('.dmnd')"
  echo
  echo "Other options:"
  echo "-b STRING         Input BAM filename"
  echo "                  Default is to run without a BAM file and therefore without expression level filtering"
  echo "-a STRING         Other argument(s) to pass to EnTAP"
  echo "-h                Print this help message and exit"
  echo
  echo "Example:    $0 -i assembly.fa -c entap_config.ini -d uniprot_sprot.dmnd -o results/entap"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "EnTAP documentation: https://entap.readthedocs.io/en/latest/"
  echo "EnTAP paper: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13106"
  echo
}

## Option defaults
assembly=""
dbs=""
outdir=""
config_file=""
bam=""
more_args=""

## Parse command-line options
while getopts 'i:d:o:c:b:a:h' flag; do
  case "${flag}" in
    i) assembly="$OPTARG" ;;
    d) dbs="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    c) config_file="$OPTARG" ;;
    b) bam="$OPTARG" ;;
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
echo "## Starting script entap.sh"
date
echo
echo "## Input assembly:                       $assembly"
echo "## Input database files:                 $dbs"
echo "## Config file:                          $config_file"
[[ "$bam" != "" ]] && echo "## Input BAM file:                       $bam"
echo "## Output dir:                           $outdir"
echo "## Other arguments to pass to EnTAP:     $more_args"
echo -e "--------------------\n"


# RUN ENTAP CONFIG -------------------------------------------------------------
echo "## Now running EnTAP..."
EnTAP \
    --runP \
    --ini "$config_file" \
    -i "$assembly" \
    $db_arg \
    --out-dir "$outdir" \
    -t "$SLURM_CPUS_PER_TASK" $more_args


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script entap.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
