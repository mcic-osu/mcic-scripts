#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --job-name=nextdenovo
#SBATCH --output=slurm-nextdenovo-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run Nextdenovo to assemble a genome."
  echo
  echo "Syntax: $0 -c <config-file>..."
  echo
  echo "Required options:"
  echo "    -c FILE           Config file for NextDenovo"
  echo "                      See the example file 'nextdenovo.cfg' in the same dir as this script"
  echo
  echo "Other options:"
  echo "    -a STRING         Other argument(s) to pass to Nextdenovo"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i data/my.fastq.gz -o results/canu -p my_genome -s 250m"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Nextdenovo documentation: https://nextdenovo.readthedocs.io/en/latest/OPTION.html"
  echo
}

## Option defaults
config_file=""
more_args=""

## Parse command-line options
while getopts ':c:a:h' flag; do
  case "${flag}" in
    c) config_file="$OPTARG" ;;
    a) more_args="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
    :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
  esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python
source activate /fs/project/PAS0471/jelmer/conda/nextdenovo-env

## Bash script settings
set -euo pipefail

## Report
echo -e "\n## Starting script nextdenovo.sh"
date
echo
echo "## Config file:                              $config_file"
[[ $more_args != "" ]] && echo "## Other arguments to pass to NextDenovo:    $more_args"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
nextDenovo \
    "$config_file" $more_args


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo -e "\n## Done with script nextdenovo.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
