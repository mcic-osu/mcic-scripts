#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=rcorrector
#SBATCH --output=slurm-rcorr-%j.out


# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "             Run Rcorrector for a directory of FASTQ files"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input dir> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  Note: Specify input either with -i/--indir or with -I/--fofn"
    echo "  -i/--indir      <file>  Input dir with FASTQ files (mutually exclusive with -I/--fofn)"
    echo "  -I/--fofn       <file>  File of file names (fofn): text file with paths to FASTQ files, one file per line"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --stage         <int>   Start at specified step           [default: '0']"
    echo "                            Should be '0' unless restarting an interrupted/failed run"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to Rcorrector"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for Rcorrector and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/fastq -o results/rcorrector"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/rcorrector-1.0.5
    set -u
}

# Print help for the focal program
Print_help_program() {
    Load_software
    run_rcorrector.pl
}

# Print SLURM job resource usage info
Resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
    echo
}

# Print SLURM job requested resources
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

# Set the number of threads/CPUs
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

# Resource usage information
Time() {
    /usr/bin/time -f \
        '\n# Ran the command:\n%C \n\n# Run stats by /usr/bin/time:\nTime: %E   CPU: %P    Max mem: %M K    Exit status: %x \n' \
        "$@"
}   

# Exit upon error with a message
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
# Option defaults
stage=0

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
indir=""
fofn=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -I | --fofn )       shift && fofn=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --more-args )       shift && more_args=$1 ;;
        --stage )           shift && stage=$1 ;;
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
# Bash script settings
set -euo pipefail

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Test parameter values
[[ "$indir" = "" && "$fofn" = "" ]] && Die "Please specify input file with -i/--indir or with -I/--fofn" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ "$indir" != "" && ! -d "$indir" ]] && Die "Input dir $indir does not exist"
[[ "$fofn" != "" && ! -f "$fofn" ]] && Die "Input fofn $fofn does not exist"

# Get list of input files
if [[ "$fofn" = "" ]]; then
    R1_list=$(echo "$indir"/*R1*fastq.gz | sed 's/ /,/g')
    R2_list=$(echo "$indir"/*R2*fastq.gz | sed 's/ /,/g')
else
    R1_list=$(grep "R1.*fastq.gz$" "$fofn" | tr "\n" ",")
    R2_list=$(grep "R2.*fastq.gz$" "$fofn" | tr "\n" ",")
fi

# Report
echo
echo "=========================================================================="
echo "                STARTING SCRIPT RCORRECTOR.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input dir:                        $indir"
echo "Output dir:                       $outdir"
[[ $more_args != "" ]] && echo "Other arguments for Rcorrector:   $more_args"
echo "Number of threads/cores:          $threads"
echo
echo "List of R1 files:                 $R1_list"
echo "List of R2 files:                 $R2_list"
echo
echo "Listing the input file(s):"
[[ $indir != "" ]] && ls -lh "$indir"
[[ $fofn != "" ]] && cat "$fofn" | xargs ls -lh
[[ $dryrun = true ]] && echo -e "\nTHIS IS A DRY-RUN"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
${e}mkdir -p "$outdir"/logs

# Run Rcorrector
${e}Time run_rcorrector.pl \
    -t "$threads" \
    -od "$outdir" \
    -stage "$stage" \
    -1 "$R1_list" \
    -2 "$R2_list" \
    $more_args


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo -e "\n# Listing files in the output dir:"
    ls -lhd "$PWD"/"$outdir"/*
    [[ "$slurm" = true ]] && Resource_usage
fi
echo "# Done with script"
date
