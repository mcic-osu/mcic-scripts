#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=inspector_correct
#SBATCH --output=slurm-inspector_correct-%j.out

#TODO - Move to output dir, or it will deposit the Bloom filter file in the working dir

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "              Run Inspector to correct a genome assembly"
    echo "======================================================================"
    echo "USAGE:"
    echo "  sbatch $0 --assembly <assembly> --reads <FASTQ> -o <output-dir> [...]"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--inspector_dir  <dir>   Dir with Inspector results"
    echo "  -o/--assembly_out   <dir>   Output dir"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --datatype          <str>   Input read type: pacbio-raw, pacbio-hifi, pacbio-corr, nano-raw, nano-corr [default: 'nano-raw']"
    echo "  --base_error                Also correct base-errors                [default: don't correct]"
    echo "  --more_args         <str>   Quoted string with additional argument(s) to pass to Inspector"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --dryrun                    Dry run: don't execute commands, only parse arguments and report"
    echo "  --debug                     Run the script in debug mode (print all code)"
    echo "  -h                          Print this help message and exit"
    echo "  --help                      Print the help for Inspector and exit"
    echo "  -v/--version                Print the version of Inspector and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --assembly results/assembly/my.fasta -o results/inspector"
    echo
    echo "SOFTWARE DOCUMENTATION:"
    echo "  - GitHub repo: https://github.com/ChongLab/Inspector / https://github.com/Maggi-Chen/Inspector"
    echo "  - Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02527-4"
    echo
}

# Load software
Load_software() {
    set +u
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    #source activate /fs/ess/PAS0471/jelmer/conda/inspector-1.2.0
    source activate /fs/ess/PAS0471/jelmer/conda/inspector-1.0.2
    set -u
}

# Print version
Print_version() {
    set +e
    Load_software
    inspector.py --version
    set -e
}

# Print help for the focal program
Print_help_program() {
    Load_software
    inspector.py --help
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
datatype="nano-raw"
base_error=false && base_error_arg="--skip_baseerror"

debug=false
dryrun=false && e=""
slurm=true


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
inspector_dir=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --inspector_dir )   shift && inspector_dir=$1 ;;
        -o | --assembly_out )    shift && assembly_out=$1 ;;
        --datatype )             shift && datatype=$1 ;;
        --more_args )            shift && more_args=$1 ;;
        --base_error )           base_error=true ;;
        -v | --version )         Print_version; exit 0 ;;
        -h )                     Print_help; exit 0 ;;
        --help )                 Print_help_program; exit 0;;
        --dryrun )               dryrun=true && e="echo ";;
        --debug )                debug=true ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
[[ -z "$SLURM_JOB_ID" ]] && slurm=false

# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Bash script settings
set -euo pipefail

# Load software and set nr of threads
[[ "$dryrun" = false ]] && Load_software
Set_threads

# Check input
[[ "$inspector_dir" = "" ]] && Die "Please specify an assembly -i/--inspector_dir"
[[ "$assembly_out" = "" ]] && Die "Please specify an output assembly file with -o/--assembly_out"
[[ ! -d "$inspector_dir" ]] && Die "Input dir (--inspector_dir) $inspector_dir does not exist"

# Determine output dir
outdir=$(dirname "$assembly_out")

# Make output assembly path absolute
[[ ! "$assembly_out" =~ ^/ ]] && assembly_out="$PWD"/"$assembly_out"

# Build other args
[[ "$base_error" = true ]] && base_error_arg=""

# Report
echo
echo "=========================================================================="
echo "          STARTING SCRIPT INSPECTOR_CORRECT.SH"
date
echo "=========================================================================="
echo "All arguments to this script:         $all_args"
echo "Input Inspector results dir:          $inspector_dir"
echo "Output assembly:                      $assembly_out"
echo "Data type:                            $datatype"
echo "Correct base-errors, too:             $base_error"
[[ $more_args != "" ]] && echo "Other arguments to pass to Inspector:     $more_args"
echo "=========================================================================="

# Print reserved resources
[[ "$slurm" = true ]] && Print_resources


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Now creating the output directories..."
${e}mkdir -pv "$outdir"/logs

# Move to the outdir or some files will go to the working dir
cd "$outdir" || exit 1

# Run
echo -e "\n# Now running Inspector-correct..."
${e}Time inspector-correct.py \
    --inspector "$inspector_dir" \
    "$base_error_arg" \
    -o . \
    --datatype "$datatype" \
    --thread "$threads" \
    $more_args

# Run
echo -e "\n# Renaming the output file..."
${e}mv -v contig_corrected.fa "$assembly_out"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
if [[ "$dryrun" = false ]]; then
    echo "# Version used:"
    Print_version | tee "$outdir"/logs/version.txt
    echo -e "\n# Listing the output assembly:"
    ls -lh "$assembly_out"
    [[ "$slurm" = true ]] && Resource_usage
fi
echo
echo "# Done with script"
date
