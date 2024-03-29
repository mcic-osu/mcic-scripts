#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=pycoqc
#SBATCH --output=slurm-pycoqc-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "    RUN PYCOQC FOR QC OF ONT DATA USING A SEQUENCING SUMMARY FILE"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> --min_pass_qual <min. qual. score> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input 'sequencing summary' file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --min_pass_qual <int>   Min. read length to PASS (should be the same as used for Guppy!)"
    echo "  --min_pass_len  <int>   Min. read length to PASS                    [default: 0]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) to pass to PycoQC"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for PycoQC and exit"
    echo "  -v/--version            Print the version of PycoQC and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i data/sequencing_summary.txt -o results/pycoqc --min_pass_qual 10"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - https://tleonardi.github.io/pycoQC/pycoQC/usage/"
    echo "  - https://tleonardi.github.io/pycoQC/pycoQC/CLI_usage/"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/pycoqc-2.5.2
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    pycoQC --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    pycoQC --help
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
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
infile=
outdir=
min_pass_qual= && qual_arg=""
min_pass_len= && len_arg=
more_args=""

# Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --min_pass_qual )   shift && min_pass_qual=$1 ;;
        --min_pass_len )    shift && min_pass_len=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        * )                 Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Bash script settings
set -euo pipefail

# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# Load software and set nr of threads
Load_software

# Check the input
[[ "$infile" = "" ]] && Die "Please specify an input file with -i/--infile" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && Die "Input file $infile does not exist"

# Other args
[[ -n "$min_pass_len" ]] && len_arg="--min_pass_len $min_pass_len"
[[ -n "$min_pass_qual" ]] && qual_arg="--min_pass_qual $min_pass_qual"

# Determine the output file name
infile_noext=${infile%.*}
outfile="$outdir"/$(basename "$infile_noext")_pycoqc.html

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT PYCOQC.SH"
date
echo "=========================================================================="
echo "All arguments to this script:     $all_args"
echo "Input file:                       $infile"
echo "Output file:                      $outfile"
[[ $more_args != "" ]] && echo "Other arguments for PycoQC:       $more_args"
echo
echo "Listing the input file(s):"
ls -lh "$infile"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
mkdir -p "$outdir"/logs

# Run
echo -e "\n# Now running PycoQC..."
Time pycoQC \
    -f "$infile" \
    -o "$outfile" \
    $qual_arg \
    $len_arg \
    $more_args


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo "# Version used:"
Print_version | tee "$outdir"/logs/version.txt
echo -e "\n# Listing the output file:"
ls -lhd "$PWD"/"$outfile"
[[ "$slurm" = true ]] && Resource_usage
echo "# Done with script"
date
