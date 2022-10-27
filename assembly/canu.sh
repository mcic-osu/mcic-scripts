#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --job-name=canu
#SBATCH --output=slurm-canu-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Print_help() {
    echo
    echo "==============================================================================="
    echo "            $0: Run Canu to assemble a genome"
    echo "==============================================================================="
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <FASTQ-file> -o <output-dir> -p <prefix> -s <genome-size> ..."
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "    -i STRING         Space-separated list of input FASTQ files"
    echo "    -o DIR            Output dir"
    echo "    -p DIR            Custom output prefix / name for the genome assembly"

    echo "    -s STRING         Estimated genome size, e.g. '65m' for 65 Mbp"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "    -t STRING         Time limit for Canu jobs: specify as HH:MM:SS     [default: 06:00:00]"
    echo "    -w DIR            Work dir - initial output dir before files are copied"
    echo "                      [default: '/fs/scratch/PAS0471/\$USER/canu/<output-prefix>']"
    echo "    -a STRING         Other argument(s) to pass to Canu"
    echo
    echo "UTILITY OPTIONS"
    echo "    -h                Print this help message and exit"
    echo "    -v                Print the version of Canu and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "    sbatch $0 -i data/my.fastq.gz -o results/canu -p my_genome -s 250m"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "- The script assumes that reads are Nanopore"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "- https://canu.readthedocs.io/en/latest/quick-start.html"
    echo
}

## Load software
Load_software() {
    module load python/3.6-conda5.2
    source activate /users/PAS0471/jelmer/miniconda3/envs/canu-env
}

## Print version
Print_version() {
    Load_software
    canu --version
}

## Exit upon error with a message
Die() {
    printf "\n$0: ERROR: %s\n" "$1" >&2
    echo -e "Exiting\n" >&2
    exit 1
}

# PARSE OPTIONS ----------------------------------------------------------------
## Constants
CNS_MEMORY=4
GRID_MEMORY="--mem=20G"

## Option defaults
WORKDIR_BASE="/fs/scratch/PAS0471/$USER/canu/"
time_limit="5:00:00"

## Placeholder defaults
infiles=""
output_prefix=""
outdir=""
genome_size=""
more_args=""
workdir=""

## Parse command-line options
while getopts ':i:o:w:p:t:a:s:vh' flag; do
    case "${flag}" in
        i) infiles="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        p) output_prefix="$OPTARG" ;;
        w) workdir="$OPTARG" ;;
        t) time_limit="$OPTARG" ;;
        s) genome_size="$OPTARG" ;;
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

## Bash script settings
set -euo pipefail

## Arguments for Canu
slurm_options="--account=$SLURM_JOB_ACCOUNT --time=$time_limit"
infile_arg="-nanopore ${infiles// / -nanopore }"

## Workdir
[[ "$workdir" = "" ]] && workdir="$WORKDIR_BASE"/"$output_prefix"

## Check input
[[ $infiles = "" ]] && Die "Please specify one or more input files with -i"
[[ $outdir = "" ]] && Die "Please specify an output dir with -o"
[[ $output_prefix = "" ]] && Die "Please specify an output prefix with -p"
[[ $genome_size = "" ]] && Die "Please specify a genome size with -s"


## Report
echo "=========================================================================="
echo "                        STARTING SCRIPT CANU.SH"
date
echo "=========================================================================="
echo "Input FASTQ file(s):                  $infiles"
echo "Output dir:                           $outdir"
echo "Output prefix:                        $output_prefix"
echo "Work dir:                             $workdir"
[[ $more_args != "" ]] && echo "Other arguments to pass to Canu:    $more_args"
echo "Listing input files:"
for infile in $infiles; do
    [[ ! -f $infile ]] && Die "Input file $infile does not exist!"
    ls -lh "$infile"
done
echo "=========================================================================="


# MAIN -------------------------------------------------------------------------
## Create the output directory
mkdir -p "$outdir"/logs

echo -e "\n## Now running Canu..."
canu \
    -p "$output_prefix" \
    -d "$workdir" \
    genomeSize="$genome_size" \
    gridOptions="$slurm_options" \
    gridEngineMemoryOption="$GRID_MEMORY" \
    cnsMemory="$CNS_MEMORY" \
    $more_args $infile_arg

##TODO - Copy some files to outdir


# WRAP UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
