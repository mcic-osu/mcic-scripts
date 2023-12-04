#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=30:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=kraken_build
#SBATCH --output=slurm-kraken_build-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Build a custom Kraken database
  Note: Standard Kraken databases can also be downloaded from https://benlangmead.github.io/aws-indexes/"
SCRIPT_VERSION="2023-12-03"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=kraken2-build
TOOL_NAME=Kraken-build
TOOL_DOCS=https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/kraken2
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                 # When true, just print tool & script version info and exit

# Defaults - tool parameters
libs=false

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
    echo "  - Basic usage:"
    echo "      sbatch $0 -o refdata/kraken/my_db --genome_dir refdata/kraken/my_db/genomes"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--db_dir         <dir>   Dir for the Kraken DB"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --libs              <str>   Comma-separated list of reference libraries to include [default:none]"
    echo "                                Options: 'archaea', 'bacteria', 'viral', 'plasmid', 'human', 'fungi', 'plants', 'protozoa',"
    echo "                                         'UniVec', 'Univec_Core', 'nr', 'nt'"
    echo "                                See https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases"
    echo "  --genome_fa         <file>  Custom genome FASTA file to be added to the db (use if wanting to add a single genome)"
    echo "  --genome_dir        <dir>   Dir with custom genome FASTA files (use if wanting to add multiple genomes)"
    echo "  --more_opts         <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                you'll have to provide one in order to run the script with a container.)"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or '--dl_container' is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: $dl_container]"
    echo "  --no_strict                 Don't use strict Bash settings ('set -euo pipefail') -- can be useful for troubleshooting"
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
    function_script=$(realpath "$script_dir"/../dev/"$(basename "$FUNCTION_SCRIPT_URL")")
    # Download the function script if needed, then source it
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        function_script=$(basename "$FUNCTION_SCRIPT_URL")
        wget -q "$FUNCTION_SCRIPT_URL" -O "$function_script"
    fi
    source "$function_script"
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
db_dir=
genome_fa=
genome_dir=
more_opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --db_dir )     shift && db_dir=$1 ;;
        --genome_fa )       shift && genome_fa=$1 ;;
        --genome_dir )      shift && genome_dir=$1 ;;
        --libs )            shift && libs=$1 ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env )             shift && env=$1 ;;
        --no_strict )       strict_bash=false ;;
        --dl_container )    dl_container=true ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 && dl_container=true ;;
        -h | --help )       script_help; exit 0 ;;
        -v )                script_version; exit 0 ;;
        --version )         version_only=true ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
[[ "$strict_bash" == true ]] && set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$db_dir" ]] && die "No db dir specified, do so with -o/--db" "$all_opts"

# Define outputs based on script parameters
LOG_DIR="$db_dir"/logs && mkdir -p "$LOG_DIR"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Output database dir:                      $db_dir"
[[ -n $genome_fa ]] && echo "Genome to be added to DB:                 $genome_fa"
[[ -n $genome_dir ]] && echo "Dir with genomes to be added to DB:       $genome_dir"
[[ $libs != false ]] && echo "Libraries to download:                    $libs"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
[[ -n $genome_fa || -n "$genome_dir" ]] && echo -e "\nListing the input genome(s):"
[[ -n $genome_fa ]] && ls -lh "$genome_fa"
[[ -n $genome_dir ]] && ls -lh "$genome_dir"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
mkdir -p "$LOG_DIR" "$db_dir"

# Download taxonomy
if [[ ! -d "$db_dir"/taxonomy ]]; then
    log_time "Downloading taxonomy..."
    runstats $TOOL_BINARY --download-taxonomy --db "$db_dir"
fi

# Download libraries
if [[ "$libs" != "false" ]]; then
    IFS="," read -ra lib_array <<< "$libs"
    for lib in "${lib_array[@]}"; do
        if [[ ! -d "$db_dir"/library/"$lib" ]]; then
            log_time "Downloading library: $lib..."
            $TOOL_BINARY --download-library "$lib" --db "$db_dir"
        else
            log_time "Library dir for $lib already exists, not downloading again"
        fi
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
log_time "Building the database..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    --build \
    --db "$db_dir" \
    -t "$threads" \
    $more_opts

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$db_dir")"/*
final_reporting "$LOG_DIR"
