#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=blastdb_download
#SBATCH --output=slurm-blastdb_download-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Download a standard NCBI BLAST database"
SCRIPT_VERSION="2025-02-16"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=update_blastdb.pl
TOOL_NAME=update_blastdb.pl
TOOL_DOCS=https://www.ncbi.nlm.nih.gov/books/NBK569850/
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env_type=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/conda/blast-2.16.0
container_dir="$HOME/containers"

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage example:"
    echo "      sbatch $0 -db nt -o results/blastdb"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --db                <file>  NCBI BLAST database to download"
    echo "                              Possible options include: nr, nt, core_nt, nt_euk, nt_viruses, swissprot, taxdb, mito"
    echo "                              Run 'update_blastdb.pl --showall' to see all options"
    echo "  -o/--outdir         <dir>   Output dir for BLAST database (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_opts         <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env_type               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env_type]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                 you'll have to provide one in order to run the script with a container.)"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
}

# Function to source the script with Bash functions
source_function_script() {
    # Determine the location of this script, and based on that, the function script
    if [[ "$IS_SLURM" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script_name="$(basename "$FUNCTION_SCRIPT_URL")"
    function_script="$script_dir"/../dev/"$function_script_name"

    # Download the function script if needed, then source it
    if [[ -f "$function_script" ]]; then
        source "$function_script"
    elif [[ ! -f "$function_script_name" ]]; then
        echo "Can't find script with Bash functions ($function_script_name), downloading from GitHub..."
        wget -q "$FUNCTION_SCRIPT_URL" -O "$function_script_name"
        source "$function_script_name"
    else
        source "$function_script_name"
    fi
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
db=
outdir=
more_opts=
threads=
container_path=
container_url=
version_only=false                 # When true, just print tool & script version info and exit

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --db )              shift && db=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env_type )             shift && env_type=$1 ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v | --version )         version_only=true ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Load software
load_env "$conda_path" "$container_path"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$db" ]] && die "No DB specified, do so with ---db" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Database to download:                     $db"
echo "Output dir:                               $outdir"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
log_time "Moving into output dir $outdir..."
cd "$outdir" || exit 1

log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --decompress \
    --num_threads "$threads" \
    $more_opts \
    "$db"

#? - To see which databases are available, run:
#?   update_blastdb.pl --showall
#? - To check which sequences are in a BLAST db, run:
#?   blastdbcmd -db <my-db> -entry all -outfmt '%a %l %T %K %S %L'

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
cd - || exit 1
log_time "Listing files in the output dir:"
ls -lh
final_reporting "$LOG_DIR"
