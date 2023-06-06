#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=kraken-build
#SBATCH --output=slurm-kraken-build-%j.out

#TODO - Process ref_libs argument!

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=kraken-build.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/kraken2-2.1.2
readonly TOOL_BINARY=kraken2-build
readonly TOOL_NAME=Kraken-build
readonly TOOL_DOCS=https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
readonly TOOL_PAPER=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0

# Option defaults
ref_libs=false

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  Build a custom Kraken database"
    echo "  Note: Standard Kraken databases can be downloaded from https://benlangmead.github.io/aws-indexes/"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-file> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--db_dir     <dir>   Dir for Kraken DB"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --ref_libs      <str>   Comma-separated list of reference libraries to include [default:none]"
    echo "                          Options: 'archaea', 'bacteria', 'viral', 'plasmid', 'human', 'fungi', 'plants', 'protozoa',"
    echo "                                    'UniVec', 'Univec_Core', 'nr', 'nt' -- see https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases"
    echo "  --genome_fa     <file>  Custom genome FASTA file to be added to the db"
    echo "  --genome_dir    <dir>   Dir with custom genome FASTA files"
    echo "  --more_args     <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for $TOOL_NAME and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -o refdata/kraken/my_db --genome_dir refdata/kraken/my_db/genomes"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
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
db_dir=""
genome_fa=""
genome_dir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --db_dir )     shift && readonly db_dir=$1 ;;
        --genome_fa )       shift && readonly genome_fa=$1 ;;
        --genome_dir )      shift && readonly genome_dir=$1 ;;
        --ref_libs )        shift && readonly ref_libs=$1 ;;
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
[[ -z "$db_dir" ]] && die "No input file specified, do so with -i/--infile" "$all_args"

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
readonly version_file="$db_dir"/logs/version.txt
readonly log_dir="$db_dir"/logs

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Output database dir:                      $db_dir"
[[ -n $genome_fa ]] && echo "Genome to be added to DB:                 $genome_fa"
[[ -n $genome_dir ]] && echo "Dir with genomes to be added to DB:       $genome_dir"
[[ -n $more_args ]] && echo "Libraries to download:                    $more_args"
[[ -n $more_args ]] && echo "Other arguments for $TOOL_NAME:           $more_args"
echo "Number of threads/cores:                  $threads"
[[ -n $genome_fa || -n "$genome_dir" ]] && echo -e "\nListing the input genome(s):"
[[ -n $genome_fa ]] && ls -lh "$genome_fa"
[[ -n $genome_dir ]] && ls -lh "$genome_dir"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir" "$db_dir"

# Download taxonomy
if [ ! -d "$db_dir"/taxonomy ]; then
    log_time "Downloading taxonomy..."
    runstats $TOOL_BINARY --download-taxonomy --db "$db_dir"
fi

# Download libraries
if [[ "$ref_libs" != "false" ]]; then
    libs=() #TODO
    for lib in "${libs[@]}"; do
        log_time "Downloading library: $lib..."
        $TOOL_BINARY --download-library "$lib" --db "$db_dir"
    done
fi

# Add custom genoms
if [[ -n "$genome_fa" ]]; then
    # If there is a single genome to be added
    log_time "Adding custom genome $genome_fa to Kraken library..."
    $TOOL_BINARY --add-to-library "$genome_fa" --db "$db_dir"
elif [[ -n "$genome_dir" ]]; then
    # If all genomes in a dir should be added
    shopt -s nullglob
    for genome_fa in "$genome_dir"/*.{fa,fasta,fna}; do
        log_time "Adding custom genome $genome_fa to Kraken library..."
        $TOOL_BINARY --add-to-library "$genome_fa" --db "$db_dir"
    done
    shopt -u nullglob
else
    log_time "Not adding any custom genomes"
fi

# Run the tool
log_time "Building the database with $TOOL_NAME..."
runstats $TOOL_BINARY \
    --build \
    --db "$db_dir" \
    -t "$threads" \
    $more_args

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
tool_version | tee "$version_file"
script_version | tee -a "$version_file" 
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$db_dir")"/*
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME"
echo

