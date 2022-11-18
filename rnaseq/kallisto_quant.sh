#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=kallisto
#SBATCH --output=slurm-kallisto-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  QUANTIFY TRANSCRIPTS WITH KALLISTO"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <R1-FASTQ-infile> -o <BAM-outdir> -r <ref-index-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  R1 FASTQ input file (The name of the R2 file will be inferred by the script.)"
    echo "  -r/--ref_index  <file>  Kallisto transcriptome index (create with 'kallisto_index.sh' script)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Kallisto"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Kallisto and exit"
    echo "  -v/--version            Print the version of Kallisto and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/trinity/assembly.fa -o results/kallisto/trans.idx"
    echo
}

## Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/kallisto-0.48.0
}

## Print version
Print_version() {
    Load_software
    kallisto version
}

## Print help for the focal program
Print_help_program() {
    Load_software
    echo "# General help:"
    kallisto -h
    echo -e "\n#Help for the 'quant' subcommand:"
    kallisto quant
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
    echo "Memory (per node):    $SLURM_MEM_PER_NODE"
    echo "CPUs per task:        $SLURM_CPUS_PER_TASK"
    [[ "$SLURM_NTASKS" != 1 ]] && echo "Nr of tasks:          $SLURM_NTASKS"
    [[ -n "$SBATCH_TIMELIMIT" ]] && echo "Time limit:           $SBATCH_TIMELIMIT"
    echo "======================================================================"
    echo
    set -u
}

## Set the number of threads/CPUs
Set_threads() {
    set +u
    if [[ "$slurm" = true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            threads="$SLURM_NTASKS"
        else 
            echo "WARNING: Can't detect nr of threads, setting to 1"
            threads=1
        fi
    else
        threads=1
    fi
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
    
    echo >&2
    echo "=====================================================================" >&2
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option" >&2
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'" >&2
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    echo -e "\nEXITING..." >&2
    echo "=====================================================================" >&2
    echo >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
## Option defaults
debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
R1=""
outdir=""
more_args=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && R1=$1 ;;
        -r | --index )      shift && index=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --more-args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --dryrun )          dryrun=true && e="echo ";;
        --debug )           debug=true ;;
        * )                 Die "Invalid option $1" "$all_args" ;;
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
Set_threads

## Bash script settings
set -euo pipefail

## Check input
[[ "$R1" = "" ]] && Die "Please specify an R1 input file with -i/--infile" "$all_args"
[[ "$index" = "" ]] && Die "Please specify an index input file with -r/--index" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$R1" ]] && Die "Input file $R1 does not exist"
[[ ! -f "$index" ]] && Die "Input file $index does not exist"

## Process parameters
R2=${R1/_R1_/_R2_}                               # R2 FASTQ file
sample_id=$(basename "$R1" | sed 's/_R1_.*//')   # Sample ID - for output files
outdir_full=$outdir/"$sample_id"


## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TODO_SCRIPTNAME"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input R1 FASTQ file:              $R1"
echo "Input R2 FASTQ file:              $R2"
echo "Output dir:                       $outdir"
[[ $more_args != "" ]] && echo "Other arguments for Kallisto:     $more_args"
echo "Sample ID (inferred by the script): $sample_id"
echo "Number of threads/cores:          $threads"
echo
echo "Listing the input file(s):"
ls -lh "$R1" "$R2" "$index"
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

## Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
${e}mkdir -p "$outdir"/logs

## Run
echo -e "\n# Now running Kallisto quant..."
${e}Time kallisto quant \
    -i "$index" \
    -o "$outdir_full" \
    -b 100 \
    "$R1" "$R2"

#? -b = number of bootstraps


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
