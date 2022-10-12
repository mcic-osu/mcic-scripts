#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --job-name=nextdenovo
#SBATCH --output=slurm-nextdenovo-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
    echo
    echo "=================================================================================================="
    echo "$0: Run Nextdenovo to assemble a genome"
    echo "=================================================================================================="
    echo
    echo "SYNTAX:"
    echo "------------------"
    echo "  sbatch $0 -c <config-file> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -c FILE           Config file for NextDenovo"
    echo "                      See the example file 'nextdenovo.cfg' in the same dir as this script."
    echo "                      This config file will set the input files, output dir, and genome size."
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -a STRING         Other argument(s) to pass to Nextdenovo"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "  sbatch 0 -c nextdenovo.cfg"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "------------------"
    echo "- https://nextdenovo.readthedocs.io/en/latest/OPTION.html"
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

## Check input
[[ $config_file = "" ]] && echo -e "\n## ERROR: Please specify a config file -c \n" >&2 && exit 1
[[ ! -f $config_file ]] && echo -e "\n## ERROR: Config file $config_file does not exist! \n" >&2 && exit 1

## Report
echo -e "\n=========================================================================="
echo "## STARTING SCRIPT NEXTDENOVO.SH"
date
echo -e "==========================================================================\n"
echo "## Config file:                              $config_file"
[[ $more_args != "" ]] && echo "## Other arguments to pass to NextDenovo:    $more_args"
echo -e "--------------------\n"

## Show config file
echo "## Printing the contents of the config file:"
cat "$config_file"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
nextDenovo "$config_file" $more_args


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo -e "\n## Done with script nextdenovo.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
