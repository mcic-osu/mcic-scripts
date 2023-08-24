#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sourmash_search
#SBATCH --output=slurm-sourmash_search-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Search a query for the most similar hit in a sourmash database"
SCRIPT_VERSION="2023-08-22"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY="sourmash"
TOOL_NAME=Sourmash
TOOL_DOCS="https://sourmash.readthedocs.io/en/latest/tutorial-basic.html#make-and-search-a-database-quickly"
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/project/PAS0471/jelmer/conda/sourmash
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                 # When true, just print tool & script version info and exit

# Defaults - tool parameters
kmer_size=31
n_results=10

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
    echo "      sbatch $0 -i results/spades/sampleA.fna --db results/sourmash_search/db/db.sbt.zip -o results/sourmash_search"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input FASTA file (the query) with .fasta / .fa / .fna extension"
    echo "  --db                <file>  Path to a Sourmash database"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --n_results         <int>   Number of results (matches) to return   [default: $n_results]"
    echo "  --kmer_size         <int>   Kmer size (should be an odd integer)    [default: $kmer_size]"
    echo "                                Will be used for creating a sketch of the query file"
    echo "                                Should match the kmer-size used for the database"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                 you'll have to provide one in order to run the script with a container.)"
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
        wget "$FUNCTION_SCRIPT_URL" -O "$function_script"
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
infile=
db=
outdir=
opts=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        --db )              shift && db=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --kmer_size )       shift && kmer_size=$1 ;;
        --n_results )       shift && n_results=$1 ;;
        --opts )            shift && opts=$1 ;;
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
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$db" ]] && die "No database file specified, do so with --db" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -f "$db" ]] && die "Database file $db does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
infile=$(realpath "$infile") # Make paths absolute because we have to move into the outdir
db=$(realpath "$db")
infile_signature=$(basename "$infile").sig
sample_id=$(basename "${infile%.*}")

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input file:                               $infile"
echo "Database:                                 $db"
echo "Kmer size:                                $kmer_size"
echo "Output dir:                               $outdir"
echo "## Sample ID:                    $sample_id"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:        $opts"
log_time "Listing the input file(s):"
ls -lh "$infile" "$db"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Move to output dir
cd "$outdir" || exit

# Create a signature for the query FASTA file
if [[ ! -f "$infile_signature" ]]; then
    log_time "Creating sourmash signature for query file $infile..."
    runstats $CONTAINER_PREFIX $TOOL_BINARY \
        sketch dna -p k="$kmer_size" "$infile"
else
    log_time "Sourmash signature for query file already exists ($infile_signature)"
fi

# Run the search
log_time "Now running sourmash search..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    search --num-results "$n_results" "$infile_signature" "$db" |
    tee > "$sample_id".txt

# Report
log_time "Showing the contents of the main output file ($PWD/$sample_id.txt):"
cat "$sample_id".txt

ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
