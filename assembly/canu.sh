#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --job-name=canu
#SBATCH --output=slurm-canu-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
    echo
    echo "=================================================================================================="
    echo "$0: Run Canu to assemble a genome"
    echo "=================================================================================================="
    echo
    echo "SYNTAX:"
    echo "------------------"
    echo "  sbatch $0 -i <FASTQ-file> -o <output-dir> -p <prefix> -s <genome-size> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "------------------"
    echo "    -i STRING         Space-separated list of input FASTQ files"
    echo "    -o DIR            Output dir"
    echo "    -p DIR            Custom output prefix / name for the genome assembly"
    echo "    -s STRING         Estimated genome size, e.g. '65m' for 65 Mbp"
    echo
    echo "OTHER OPTIONS:"
    echo "------------------"
    echo "    -a STRING         Time limit for Canu jobs, specify as HH:MM:SS     [default: 06:00:00]"
    echo "    -a STRING         Other argument(s) to pass to Canu"
    echo "    -h                Print this help message and exit"
    echo
    echo "EXAMPLE COMMAND:"
    echo "------------------"
    echo "    sbatch $0 -i data/my.fastq.gz -o results/canu -p my_genome -s 250m"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "------------------"
    echo "- https://canu.readthedocs.io/en/latest/quick-start.html"
    echo
}

## Option defaults
time_limit="06:00:00"
infiles=""
output_prefix=""
outdir=""
genome_size=""
more_args=""

## Parse command-line options
while getopts ':i:o:p:t:a:s:h' flag; do
    case "${flag}" in
        i) infiles="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        p) output_prefix="$OPTARG" ;;
        t) time_limit="$OPTARG" ;;
        s) genome_size="$OPTARG" ;;
        a) more_args="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/canu-env

## Bash script settings
set -euo pipefail

## Constants
CNS_MEMORY=4
GRID_MEMORY="--mem=20G"

## Arguments for Canu
slurm_options="--account=$SLURM_JOB_ACCOUNT --time=$time_limit"
infile_arg="-nanopore ${infiles// / -nanopore }"

## Check input
[[ $infiles = "" ]] && echo -e "\n## ERROR: Please specify one or more input files with -i \n" >&2 && exit 1
[[ $outdir = "" ]] && echo -e "\n## ERROR: Please specify an output dir with -o \n" >&2 && exit 1
[[ $output_prefix = "" ]] && echo -e "\n## ERROR: Please specify an output prefix with -p \n" >&2 && exit 1
[[ $genome_size = "" ]] && echo -e "\n## ERROR: Please specify a genome size with -s \n" >&2 && exit 1


## Report
echo -e "\n=========================================================================="
echo "## STARTING SCRIPT CANU.SH"
date
echo -e "==========================================================================\n"
echo "## Input FASTQ file(s):                  $infiles"
echo "## Output dir:                           $outdir"
echo "## Output prefix:                        $output_prefix"
[[ $more_args != "" ]] && echo "## Other arguments to pass to Canu:    $more_args"
echo "## Listing input files:"
for infile in $infiles; do
    [[ ! -f $infile ]] && echo -e "\n## ERROR: Input file $infile does not exist! \n" >&2 && exit 1
    ls -lh "$infile"
done
echo -e "--------------------\n"


# MAIN -------------------------------------------------------------------------
## Create the output directory
mkdir -p "$outdir"

echo "## Now running Canu..."
canu \
    -p "$output_prefix" \
    -d "$outdir" \
    genomeSize="$genome_size" \
    gridOptions="$slurm_options" \
    gridEngineMemoryOption="$GRID_MEMORY" \
    cnsMemory="$CNS_MEMORY" \
    $infile_arg $more_args


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script canu.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
