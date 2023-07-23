#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=panaroo
#SBATCH --output=slurm-panaroo-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run a pangenome analysis with Panaroo"
MODULE=miniconda3
#CONDA=/fs/ess/PAS0471/jelmer/conda/panaroo #! Getting an error with this!
CONTAINER=/fs/ess/PAS0471/containers/depot.galaxyproject.org-singularity-panaroo-1.3.3--pyhdfd78af_0.img
SCRIPT_VERSION="2023-07-22"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY="singularity exec $CONTAINER panaroo"
TOOL_NAME=Panaroo
TOOL_DOCS=https://gtonkinhill.github.io/panaroo
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - settings
ALIGN_THESE=core            # Align 'core' or 'pan' (all genes); Panaroo's '--alignment' argument

# Parameter defaults
mode=strict                 # Use Panaroo's 'strict' mode by default
core_threshold=0.95         # Same as Panaroo's default

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
    echo "      sbatch $0 -i results/prokka -o results/panaroo"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir      <file>  Input dir with Prokka or Bakta GFF files ('.gff' or '.gff3' extension)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --clean_mode    <str>   Panaroo's running mode: 'strict', 'moderate', or 'sensitive' [default: 'strict']"
    echo "  --core_threshold <num>  Core-genome sample threshold                [default=0.95]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "HARDCODED OPTIONS:"
    echo "  The script always uses '--alignment core' to only align core genes"
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
indir=
outdir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --clean_mode )      shift && mode=$1 ;;
        --core_threshold )  shift && core_threshold=$1 ;;
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
[[ -z "$indir" ]] && die "No input file specified, do so with -i/--indir" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -d "$indir" ]] && die "Input dir $indir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Logging files and dirs
LOG_DIR="$outdir"/logs
VERSION_FILE="$LOG_DIR"/version.txt
#CONDA_YML="$LOG_DIR"/conda_env.yml
ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
#load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input dir:                                    $indir"
echo "Number of GFF files:                          $(ls "$indir"/*gff? | wc -l)"
echo "Output dir:                                   $outdir"
echo "Panaroo mode:                                 $mode"
echo "Align 'pan' or 'core' genes:                  $ALIGN_THESE" 
echo "Core gene threshold:                          $core_threshold"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$indir"/*gff?
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    -i "$indir"/*gff? \
    -o "$outdir" \
    --clean-mode "$mode" \
    --core_threshold "$core_threshold" \
    --alignment "$ALIGN_THESE" \
    --threads "$threads" \
    $more_args

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
