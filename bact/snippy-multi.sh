#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=snippy
#SBATCH --output=slurm-snippy-%j.out

# FUNCTIONS --------------------------------------------------------------------
# Help function
print_help() {
    echo
    echo "======================================================================"
    echo                         "$0"
    echo     "Run snippy-multi to align FASTQ files for multiple samples"
    echo "                to a reference genome and find SNPs"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-table> -r <reference> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input TSV with paths to FASTQ files"
    echo "                          Should have 3 columns (no header): genome ID, forward reads FASTQ, reverse reads FASTQ"
    echo "  -r/--ref        <file>  Input reference genome FASTA or Genbank file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to snippy-multi"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for snippy-multi and exit"
    echo "  -v/--version            Print the version of Snippy and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i snippy-input.tsv -r results/assembly/genome.fa -o results/snippy"
    echo
    echo "DOCUMENTATION:"
    echo "  - Repo/documentation:   https://github.com/tseemann/snippy"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/snippy-4.6.0
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    snippy --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    snippy-multi --help
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
    echo "Memory (MB per node): $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):      $SLURM_CPUS_PER_TASK"
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
# Option defaults
debug=false
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
input_table=""
ref=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --input_table )    shift && input_table=$1 ;;
        -r | --ref )            shift && ref=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        -v | --version )        Print_version; exit 0 ;;
        -h )                    Print_help; exit 0 ;;
        --help )                Print_help_program; exit 0;;
        --debug )               debug=true ;;
        * )                     Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
Load_software
Set_threads

# Check input
[[ "$input_table" = "" ]] && Die "Please specify an input table with -i/--input_table" "$all_args"
[[ "$ref" = "" ]] && Die "Please specify an reference genome file -r/--ref" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$input_table" ]] && Die "Input file $input_table does not exist"
[[ ! -f "$ref" ]] && Die "Input file $ref does not exist"

# Make path to reference absolute
[[ ! "$ref" =~ ^/ ]] && ref="$PWD"/"$ref"

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT SNIPPY-MULTI.SH"
date
echo "=========================================================================="
echo "All arguments to this script:         $all_args"
echo "Input TSV with paths to FASTQ files:  $input_table"
echo "Reference genome file:                $ref"
echo "Output dir:                           $outdir"
[[ $more_args != "" ]] && echo "Other arguments for Snippy:           $more_args"
echo "Number of threads/cores:              $threads"
echo
echo "# Listing the input file(s):"
ls -lh "$input_table" "$ref"
echo
echo "# Showing the contents of the input table file:"
cat -n "$input_table"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
mkdir -pv "$outdir"/logs

echo -e "\n# Generating Snippy commands with snippy-multi..."
Time snippy-multi \
    "$input_table" \
    --ref "$ref" \
    --cpus "$threads" \
    --force \
    $more_args \
    > "$outdir"/runme.sh

echo -e "\n# Showing contents of the generated runme.sh file:"
cat -n "$outdir"/runme.sh

echo -e "\n# Now running Snippy..."
cd "$outdir" || exit 1
Time bash runme.sh


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing files in the output dir:"
ls -lhd "$PWD"/*
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
echo
