#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=dl-fq
#SBATCH --output=slurm-dl-fq-%j.out

# Download FASTQ files from SRA/ENA

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=dl-SRA-fq.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/fastq-dl
readonly TOOL_BINARY=fastq-dl
readonly TOOL_NAME=fastq-dl
readonly TOOL_DOCS=https://github.com/rpetit3/fastq-dl

# Option defaults
unzip=false

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  Download FASTQ files from SRA/ENA with fastq-dl"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  sbatch $0 -a SRR5506722 -o data/sra"
    echo "  sbatch $0 -a SRR5506722,SRR6942483 -o data/sra"
    echo "  sbatch $0 -a data/sra/accessions.txt -o data/sra"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -a/--accessions <str>   Comma-separated list of one or more SRA accession numbers,"
    echo "                          or a file with accession numbers, one per line."
    echo "  -o/--outdir      <dir>  Output directory"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --unzip                 Unzip the downloaded FASTQ files            [default: keep gzipped]"
    echo "  --more_args     <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for $TOOL_NAME and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "OUTPUT:"
    echo "  - One or more, optionally gzipped, FASTQ files in the specified output directory" 
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo
}

# Load software
load_tool_conda() {
    set +u
    module load "$MODULE" # Load the OSC Conda module
    # Deactivate any active Conda environments:
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi
    source activate "$CONDA_ENV" # Activate the focal environment
    set -u
}

# Exit upon error with a message
die() {
    local error_message=${1}
    local error_args=${2-none}
    log_time "$0: ERROR: $error_message" >&2
    log_time "For help, run this script with the '-h' option" >&2
    if [[ "$error_args" != "none" ]]; then
        log_time "Arguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    log_time "EXITING..." >&2
    exit 1
}

# Log messages that include the time
log_time() { echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')]" ${1-""}; }

# Print the script version
script_version() {
    echo "Run using $SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION ($SCRIPT_URL)"
}

# Print the tool's version
tool_version() {
    set +e
    load_tool_conda
    "$TOOL_BINARY" --version
    set -e
}

# Print the tool's help
tool_help() {
    load_tool_conda
    "$TOOL_BINARY" --help
}

# Print SLURM job resource usage info
resource_usage() {
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
}

# Print SLURM job requested resources
slurm_resources() {
    set +u
    log_time "SLURM job information:"
    echo "Account (project):                        $SLURM_JOB_ACCOUNT"
    echo "Job ID:                                   $SLURM_JOB_ID"
    echo "Job name:                                 $SLURM_JOB_NAME"
    echo "Memory (MB per node):                     $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):                          $SLURM_CPUS_PER_TASK"
    echo "Time limit:                               $SLURM_TIMELIMIT"
    echo -e "=================================================================\n"
    set -u
}

# Set the number of threads/CPUs
set_threads() {
    set +u
    if [[ "$is_slurm" == true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            readonly threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            readonly threads="$SLURM_NTASKS"
        else 
            log_time "WARNING: Can't detect nr of threads, setting to 1"
            readonly threads=1
        fi
    else
        readonly threads=1
    fi
    set -u
}

# Resource usage information for any process
runstats() {
    /usr/bin/time -f \
        "\n# Ran the command: \n%C
        \n# Run stats by /usr/bin/time:
        Time: %E   CPU: %P    Max mem: %M K    Exit status: %x \n" \
        "$@"
}

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
accessions=
outdir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -a | --accessions ) shift && readonly accessions=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        --unzip )           shift && unzip=true ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -v )                script_version; exit 0 ;;
        -h )                script_help; exit 0 ;;
        --version )         tool_version; exit 0 ;;
        --help )            tool_help; exit 0;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$accessions" ]] && die "No accessions specified, do so with -a/--accessions" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
load_tool_conda
set_threads

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Define outputs based on script parameters
readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs

# Getting the accessions
if [[ ! -f "$accessions" ]]; then
    IFS=',' read -ra accession_array <<< "$accessions"
else
    mapfile -t accession_array <"$accessions"
fi

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Output dir:                               $outdir"
[[ -n $more_args ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
echo "Number of threads/cores:                  $threads"
[[ -f "$accessions" ]] && echo "Accessions file:                          $accessions"
echo "Number of accessions:                     ${#accession_array[@]}"
echo "List of accessions:                       ${accession_array[*]}"
echo
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir"

# Run the tool
log_time "Starting downloads..."
for accession in "${accession_array[@]}"; do
    log_time "Now downloading accession $accession"
    runstats $TOOL_BINARY \
        --accession "$accession" \
        --outdir "$outdir" \
        --cpus "$threads" \
        $more_args
done

if [[ "$unzip" == true ]]; then
    log_time "Unzipping FASTQ files..."
    gunzip -v "$outdir"/*gz
fi

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
tool_version | tee "$version_file"
script_version | tee -a "$version_file" 
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME"
echo

# Alternative: Use sra-tools
#/fs/ess/PAS0471/jelmer/conda/sra-tools
#prefetch "$SRA_ID" -O "$outdir"
#fasterq-dump "$SRA_ID" -O "$outdir"
