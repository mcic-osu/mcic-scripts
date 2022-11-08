#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=12
#SBATCH --time=36:00:00
#SBATCH --job-name=smartdenovo
#SBATCH --output=slurm-smartdenovo-%j.out

# FUNCTIONS ---------------------------------------------------------------------
## Help function
function Print_help() {
    echo
    echo "==========================================================================="
    echo "              $0: Run SmartDenovo to assemble a genome"
    echo "==========================================================================="
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <FASTQ-file(s)> -o <output-prefix> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i FILE/STRING    One or more input FASTQ files, quoted and space-separated"
    echo "                      For example: '-i \"data/fastq/A.fq.gz data/fastq/B.fq.gz\"'"
    echo "    -o STRING         Output prefix: output dir + genome ID"
    echo "                      For example: '-o results/smartdenovo/sus_scrofa'"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -l INTEGER        Minimum read length, shorter reads will be removed     [default: 5000]"
    echo "    -a STRING         Other argument(s) to pass to SmartDenovo"
    echo
    echo "UTILITY OPTIONS:"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i \"data/fastq/A.fq.gz data/fastq/B.fq.gz\" -o results/smartdenovo/sus_scrofa"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "- https://github.com/ruanjue/smartdenovo "
    echo
}

## Load software
Load_software() {
    module load python/3.6-conda5.2
    source activate /fs/project/PAS0471/jelmer/conda/smartdenovo-env
}

## Print version
Print_version() {
    Load_software
    smartdenovo.pl --version
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}

# PARSE OPTIONS ----------------------------------------------------------------
## Option defaults
min_readlen=5000
infiles=""
output_prefix=""
more_args=""

## Parse command-line options
while getopts 'i:o:l:a:vh' flag; do
    case "${flag}" in
        i) infiles="$OPTARG" ;;
        o) output_prefix="$OPTARG" ;;
        l) min_readlen="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        v) Print_version; exit 0 ;;
        h) Print_help; exit 0 ;;
        \?) Print_help; Die "Invalid option $OPTARG" ;;
        :) Print_help; Die "Option -$OPTARG requires an argument" ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
Load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ $infiles = "" ]] && Die "Please specify one or more input files with -i"
[[ $output_prefix = "" ]] && Die "Please specify an output prefix with -o"

## Report
echo "=========================================================================="
echo "                       STARTING SCRIPT SMARTDENOVO.SH"
date
echo "=========================================================================="
echo "Input files:           $infiles"
echo "Output prefix:         $output_prefix"
echo "Min. read length:      $min_readlen"
[[ $more_args != "" ]] && echo "Other arguments to pass to Smartdenovo:    $more_args"
echo
echo "Listing input files:"
for infile in $infiles; do
    [[ ! -f $infile ]] && Die "Input file $infile does not exist!"
    ls -lh "$infile"
done
echo "=========================================================================="


# MAIN -------------------------------------------------------------------------
## Create output dir
outdir=$(dirname "$output_prefix")
mkdir -p "$outdir"

## Generate a Makefile for smartdenovo to run
smartdenovo.pl \
    -p "$output_prefix" \
    -t "$SLURM_CPUS_PER_TASK" \
    -c 1 \
    -J "$min_readlen" \
    $infiles \
    > "$output_prefix".mak

#? "After assembly, the raw unitigs are reported in file prefix.lay.utg and consensus unitigs in prefix.cns"
#? -c 1 => make consensus
#? -J min read length -- default = 5,000

## Run SmartDenovo by running the Makefile
make -f "$output_prefix".mak


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo "## Listing files in the output dir:"
ls -lhd "$PWD"/"$outdir"/*
echo -e "\n## Done with script"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
