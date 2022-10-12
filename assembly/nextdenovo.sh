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
    echo "    -i FILE           File of file names (fofn) of input FASTQ files"
    echo "    -o DIR            Output directory"
    echo "    -s STRING         Estimated genome size"
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
fofn=""
outdir=""
genome_size=""
more_args=""

## Parse command-line options
while getopts ':i:o:s:a:h' flag; do
    case "${flag}" in
        i) fofn="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
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
source activate /fs/project/PAS0471/jelmer/conda/nextdenovo-env
MCIC_SCRIPTS_REPO=https://github.com/mcic-osu/mcic-scripts.git
TEMPLATE_CONFIG=mcic-scripts/assembly/nextdenovo.cfg

## Bash script settings
set -euo pipefail

## Check input
[[ $fofn = "" ]] && echo -e "\n## ERROR: Please specify a file of file names with -i \n" >&2 && exit 1
[[ $outdir = "" ]] && echo -e "\n## ERROR: Please specify an output dir -o \n" >&2 && exit 1
[[ $genome_size = "" ]] && echo -e "\n## ERROR: Please specify an estimated genome size with -s \n" >&2 && exit 1
[[ ! -f $fofn ]] && echo -e "\n## ERROR: File $fofn does not exist! \n" >&2 && exit 1

## Make path to FOFN absolute, if needed
[[ ! "$fofn" =~ ^/ ]] && fofn="$PWD"/"$fofn"


## Report
echo -e "\n=========================================================================="
echo "## STARTING SCRIPT NEXTDENOVO.SH"
date
echo -e "==========================================================================\n"
echo "## File of file names:                       $fofn"
echo "## Output dir:                               $outdir"
echo "## Estimated genome size:                    $genome_size"
[[ $more_args != "" ]] && echo "## Other arguments to pass to NextDenovo:    $more_args"
echo -e "--------------------\n"
echo "## Printing the contents of the fofn file:"
echo -e "--------------------"
cat -n "$fofn"
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
## Create config file
echo "## Creating the config file for Nextdenovo..."
config_file="$outdir"/nextdenovo.cfg
[[ ! -f "$TEMPLATE_CONFIG" ]] && git clone "$MCIC_SCRIPTS_REPO"
cp -v "$TEMPLATE_CONFIG" "$config_file"

sed -i  -e "s@/absolute/path/to/file.fofn@$fofn@" \
        -e "s@path/to/workdir@$outdir@" \
        -e "s/XXm/$genome_size/" \
        "$config_file"

## Show config file contents
echo "## Printing the contents of the config file:"
echo -e "--------------------"
cat -n "$config_file"
echo -e "--------------------\n"

## Make the output dir
mkdir -p "$outdir"

## Run nextDenovo
nextDenovo "$config_file" $more_args


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo -e "\n## Done with script nextdenovo.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
