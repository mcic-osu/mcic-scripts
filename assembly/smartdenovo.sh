#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --job-name=smartdenovo
#SBATCH --output=slurm-smartdenovo-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
    echo
    echo "=================================================================================================="
    echo "$0: Run SmartDenovo to assemble a genome"
    echo "=================================================================================================="
    echo
    echo "SYNTAX:"
    echo "------------------"
    echo "  sbatch $0 -i <FASTQ-file(s)> -o <output-prefix> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i FILE/STRING    One or more input FASTQ files, quoted and space-separated"
    echo "                      For example: '-i \"data/fastq/A.fq.gz data/fastq/B.fq.gz\"'"
    echo "    -o STRING         Output prefix: output dir + genome ID"
    echo "                      For example: '-o results/smartdenovo/sus_scrofa'"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -a STRING         Other argument(s) to pass to SmartDenovo"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "  sbatch $0 -i \"data/fastq/A.fq.gz data/fastq/B.fq.gz\" -o results/smartdenovo/sus_scrofa"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "------------------"
    echo "- https://github.com/ruanjue/smartdenovo "
    echo
}

## Option defaults
infiles=""
output_prefix=""
more_args=""

## Parse command-line options
while getopts 'i:o:a:h' flag; do
    case "${flag}" in
        i) infiles="$OPTARG" ;;
        o) output_prefix="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/smartdenovo-env

## Bash strict settings
set -euo pipefail

## Check input
[[ $infiles = "" ]] && echo -e "\n## ERROR: Please specify one or more input files with -i \n" >&2 && exit 1
[[ $output_prefix = "" ]] && echo -e "\n## ERROR: Please specify an output prefix with -o \n" >&2 && exit 1

## Report
echo -e "\n=========================================================================="
echo "## STARTING SCRIPT SMARTDENOVO.SH"
date
echo -e "==========================================================================\n"
echo "## Input files:       $infiles"
echo "## Output prefix:     $output_prefix"
[[ $more_args != "" ]] && echo "## Other arguments to pass to Smartdenovo:    $more_args"
echo
echo "## Listing input files:"
for infile in $infiles; do
    [[ ! -f $infile ]] && echo -e "\n## ERROR: Input file $infile does not exist! \n" >&2 && exit 1
    ls -lh "$infile"
done
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
## Create output dir
outdir=$(dirname "$output_prefix")
mkdir -p "$outdir"

## Generate a Makefile for smartdenovo to run
smartdenovo.pl \
    -p "$output_prefix" \
    -t "$SLURM_CPUS_PER_TASK" \
    -c 1 \
    $infiles \
    > "$output_prefix".mak

#? "After assembly, the raw unitigs are reported in file prefix.lay.utg and consensus unitigs in prefix.cns"
#? -c 1 => make consensus
#! -J min read length -- default = 5,000

## Run SmartDenovo by running the Makefile
make -f "$output_prefix".mak


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Done with script smartdenovo.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
