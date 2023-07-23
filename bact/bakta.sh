#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=bakta
#SBATCH --output=slurm-bakta-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Bakta to annotate a bacterial genome"
MODULE=miniconda3
CONDA=/fs/ess/PAS0471/jelmer/conda/bakta
SCRIPT_VERSION="2023-07-22"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY=bakta
TOOL_NAME=Bakta
TOOL_DOCS=https://github.com/oschwengers/bakta
VERSION_COMMAND="$TOOL_BINARY --version"

# Constants - parameters
DB_TYPE=full                    # Full rather than partial Bakta DB

# Parameter defaults
db_dir=/fs/ess/PAS0471/jelmer/refdata/bakta
force_db_download=false

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage:"
    echo "      sbatch $0 -i results/spades/sampleA.fna -o results/bakta --genus Salmonella --species enterica"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input assembly FASTA file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "  --genus         <str>   Genus name of focal organism"
    echo "  --species       <str>   Species name of focal organism"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --db_dir        <dir>   Dir with/for the Bakta DB (script will download the DB if dir doesn't exist)"
    echo "  --force_db_dl           Force Bakta DB download, even if DB dir exists"
    echo "  --min_contig_len <int>  Minimum contig length                       [default: Bakta default, 1 as of writing]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
    echo
}

# Function to source the script with Bash functions
source_function_script() {
    local is_slurm=$1

    # Determine the location of this script, and based on that, the function script
    if [[ "$is_slurm" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script=$(realpath "$script_dir"/../dev/bash_functions.sh)
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions.sh
    fi
    source "$function_script"
}

# ==============================================================================
#                          INFRASTRUCTURE SETUP I
# ==============================================================================
# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script $IS_SLURM

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
infile=
outdir=
genus=
species=
min_contig_len= && contig_len_arg=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --genus )           shift && genus=$1 ;;
        --species )           shift && species=$1 ;;
        --db_dir )          shift && db_dir=$1 ;;
        --force_db_dl )     force_db_download=true ;;
        --min_contig_len )  min_contig_len=$1 ;;
        --more_args )       shift && more_args=$1 ;;
        -v )                script_version; exit 0 ;;
        -h | --help )       script_help; exit 0 ;;
        --version )         load_env "$MODULE" "$CONDA"
                            tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$genus" ]] && die "No genus name specified, do so with --genus" "$all_args"
[[ -z "$species" ]] && die "No species name specified, do so with --species" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Logging files and dirs
LOG_DIR="$outdir"/logs
VERSION_FILE="$LOG_DIR"/version.txt
CONDA_YML="$LOG_DIR"/conda_env.yml
ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# Define outputs and settings based on provided options
[[ -n "$min_contig_len" ]] && contig_len_arg="--min-contig-length $min_contig_len"
infile_base=$(basename "$infile")
sample_id=${infile_base%.*}

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input FASTA file:                             $infile"
echo "Output dir:                                   $outdir"
echo "Genus:                                        $genus"
echo "Species:                                      $species"
echo "Bakta DB dir:                                 $db_dir"
echo "Force DB download:                            $force_db_download"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Download the Bakta DB
if [[ ! -d "$db_dir"/db || "$force_db_download" == true ]]; then
    log_time "Downloading the Bakta DB to $db_dir"
    bakta_db download --output "$db_dir" --type "$DB_TYPE"
else
    log_time "Found a Bakta DB in $db_dir/db, using that one..."
fi

# Run Bakta
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --threads "$threads" \
    --output "$outdir" \
    --prefix "$sample_id" \
    --genus "$genus" \
    --species "$species" \
    --db "$db_dir"/db \
    $contig_len_arg \
    --force \
    $more_args \
    "$infile"

#? Other options
#--compliant           Force Genbank/ENA/DDJB compliance
#  --translation-table {11,4}
#                        Translation table: 11/4 (default = 11)
#  --gram {+,-,?}        Gram type for signal peptide predictions: +/-/? (default = ?)

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
