#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
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
    echo "  (Choose one of the input file options: -i or -I)"
    echo "  -i/--infiles    <file>  An input FASTQ file"
    echo "                            Or optionally multiple files, quoted and space-separated"
    echo "  -I/--fofn       <file>  Text file with list of input FASTQ files one per line (fofn)"
    echo "  -o/--outdir     <dir>   Output dir - a scratch dir is recommended,"
    echo "                             since there can be a lot of output, most of which is not useful."
    echo "  -p/--out_prefix <str>   Prefix for output files (e.g. genome ID)"
    echo "  --genome_size   <str>   Genome size estimate, e.g '4.6m' or '1g'"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --mhap_sensitivity <str> 'low', 'normal', or 'high'                 [default: 'normal']"
    echo "                          Based on read coverage:"
    echo "                            - 'low' sensitivity is used if coverage is more than 60"
    echo "                            - 'normal' is used if coverage is between 60 and 30"
    echo "                            - 'high' is used for coverages less than 30."
    echo "  --fast                  Turn on 'fast' mode, should work for genomes <1GB [default: off]"
    echo "  --time          <str>   Time limit for Canu jobs: specify as HH:MM:SS     [default: 12:00:00]"
    echo "  --more_args     <str>   Other argument(s) to pass to Canu as a quoted string"
    echo
    echo "UTILITY OPTIONS"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Canu and exit"
    echo "  -v/--version            Print the version of Canu and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/my.fastq.gz -o /fs/scratch/PAS0471/$USER/canu --genome_size 250m"
    echo "  sbatch $0 -i \"data/fastq/A.fq.gz data/fastq/B.fq.gz\" -o /fs/scratch/PAS0471/$USER/canu --genome_size 1g"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - The script assumes that reads are Nanopore"
    echo "  - The script will have Canu submit jobs to the SLURm queue"
    echo
    echo "NOTES"
    echo "  - This script will only run for a few minutes: Canu will submit its own SLURM jobs"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - https://canu.readthedocs.io/en/latest/quick-start.html"
    echo
}

## Load software
Load_software() {
    module load python/3.6-conda5.2
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /users/PAS0471/jelmer/miniconda3/envs/canu-env
}

## Print version
Print_version() {
    Load_software
    canu -version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    canu --help
}

## Print SLURM job resource usage info
Resource_usage() {
    echo
    ${e}sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime,MaxVMSize | \
        grep -Ev "ba|ex"
    echo
}

## Print SLURM job requested resources
Print_resources() {
    set +u
    echo "# SLURM job information:"
    echo "Account (project):    $SLURM_JOB_ACCOUNT"
    echo "Job ID:               $SLURM_JOB_ID"
    echo "Job name:             $SLURM_JOB_NAME"
    echo "Memory (MB per node): $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

## Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

## Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo
    echo "====================================================================="
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option"
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'"
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:"
        echo "$error_args"
    fi
    echo -e "\nEXITING..." >&2
    echo "====================================================================="
    echo
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Constants
SLURM_ARG_INIT="--account=$SLURM_JOB_ACCOUNT"

## Option defaults
mhap_sensitivity="normal"
fast=false && fast_arg=""
time_limit="12:00:00"

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
declare -a infiles
fofn=""
outdir=""
genome_size=""
more_args=""
infile_arg=""

## Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infiles )        shift && infiles=($1) ;;
        -I | --fofn )           shift && fofn=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        -p | --out_prefix )     shift && out_prefix=$1 ;;
        --time )                shift && time_limit=$1 ;;
        --genome_size )         shift && genome_size=$1 ;;
        --mhap_sensitivity )    shift && mhap_sensitivity=$1 ;;
        --fast )                fast=true ;;
        --more_args )           shift && more_args=$1 ;;
        -v | --version )        Print_version; exit 0 ;;
        -h )                    Print_help; exit 0 ;;
        --help )                Print_help_program; exit 0 ;;
        --dryrun )              dryrun=true && e="echo " ;;
        --debug )               debug=true ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

## Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software

## Bash strict settings
set -euo pipefail

## Check input
[[ ${#infiles[@]} = 0 && "$fofn" = "" ]] && Die "Please specify input files with -i or -I" "$all_args"
[[ $outdir = "" ]] && Die "Please specify an output file with -o" "$all_args"
[[ $out_prefix = "" ]] && Die "Please specify an output prefix with -p" "$all_args"
[[ $genome_size = "" ]] && Die "Please specify a genome size with --genome_size" "$all_args"

## If a FOFN was provided, read file list into an array
[[ "$fofn" != "" ]] && mapfile -t infiles <"$fofn"

## Build the input file arg
for infile in "${infiles[@]}"; do
    infile_arg="$infile_arg -nanopore $infile"
done

## Full SLURM options
slurm_arg="$SLURM_ARG_INIT --time=$time_limit"

## Fast arg
[[ "$fast" = true ]] && fast_arg="--fast"

## Report
echo "=========================================================================="
echo "                      STARTING SCRIPT CANU.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo
[[ "$fofn" != "" ]] && echo "File with list of FASTQs (fofn):      $fofn"
echo "Input files:                          ${infiles[*]}"
echo "Output dir:                           $outdir"
echo "Output prefix:                        $out_prefix"
echo "Genome size:                          $genome_size"
echo "Time limit for Canu jobs:             $time_limit"
echo "Fast mode:                            $fast"
echo "Mhap sensitivity:                     $mhap_sensitivity"
echo
echo "Number of input files:                ${#infiles[*]}"
echo "Input file argument:                  $infile_arg"
echo "Output file prefix:                   $out_prefix"
echo "Full SLURM options for Canu jobs:     $slurm_arg"
[[ $more_args != "" ]] && echo "Other arguments for Canu:      $more_args"
echo
echo "Listing input files:"
for infile in "${infiles[@]}"; do
    [[ ! -f $infile ]] && Die "Input file $infile does not exist!"
    ls -lh "$infile"
done
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
if [[ "$dryrun" = false ]]; then

    ## Create the output directory
    mkdir -p "$outdir"/logs

    echo "## Now running Canu..."
    Time canu \
        -p "$out_prefix" \
        -d "$outdir" \
        genomeSize="$genome_size" \
        MhapSensitivity="$mhap_sensitivity" \
        gridEngine=slurm \
        gridOptions="$slurm_arg" \
        executiveMemory=8 \
        ovsMemory=40 \
        merylMemory=40 \
        merylThreads=10 \
        minMemory=40 \
        minThreads=10 \
        maxMemory=175 \
        maxThreads=42 \
        stageDirectory=\$TMPDIR \
        $more_args $fast_arg $infile_arg

fi

## Options used previously:
#CNS_MEMORY=4 
#GRID_MEMORY="--mem=40G"
#cnsMemory="$CNS_MEMORY" gridEngineMemoryOption="$GRID_MEMORY"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
