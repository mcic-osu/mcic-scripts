#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=treetime
#SBATCH --output=slurm-treetime-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run TreeTime to date a phylogenetic tree"
MODULE=miniconda3
CONDA=/fs/project/PAS0471/jelmer/conda/treetime
SCRIPT_VERSION="2023-07-21"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY=treetime
TOOL_NAME=Treetime
TOOL_DOCS=https://treetime.readthedocs.io/
VERSION_COMMAND="$TOOL_BINARY --version"

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
    echo "      sbatch $0 -i results/snippy/core.full.aln -o results/treetime --dates metadata/dates.tsv"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--alignment  <file>  Input alignment file in FASTA, Phylip, or VCF format"
    echo "  --dates         <file>  Input CSV/TSV with sampling date for each sample"
    echo "                            Should have columns 'node_name' and 'date',"
    echo "                            with dates in %Y-%m-%d format"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --tree          <file>  Tree file in Nexus or Newick format"
    echo "                            If no tree file is provided, TreeTime will infer one"
    echo "  --clock_rate    <str>   Clock rate to assume"
    echo "  --more_args     <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - $TOOL_DOCS"
    echo "  - Repo:                 https://github.com/neherlab/treetime"
    echo "  - Paper:                https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5758920/"
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
date_file=
tree= && tree_arg=
clock_rate= && clock_rate_arg=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --dates )           shift && date_file=$1 ;;
        --tree )            shift && tree=$1 ;;
        --clock_rate )      shift && clock_rate=$1 ;;
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
[[ -z "$infile" ]] && die "No input alignment file specified, do so with -i/--infile" "$all_args"
[[ -z "$date_file" ]] && die "No input dates file specified, do so with -i/--dates" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -f "$date_file" ]] && die "Input file $date_file does not exist"
[[ -n "$tree" && ! -f "$tree" ]] && die "Input tree file $tree does not exist"

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

# Define outputs based on script parameters
[[ -n $tree ]] && tree_arg="--tree $tree"
[[ -n $clock_rate ]] && clock_rate_arg="--clock-rate $clock_rate"

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input file:                                   $infile"
echo "Output dir:                                   $outdir"
echo "CSV/TSV dates input file:                     $date_file"
[[ -n "$tree" ]] && echo "Input tree file:                              $tree"
[[ -n "$clock_rate" ]] && echo "Clock rate:                                   $clock_rate"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$infile" "$date_file"
[[ -n "$tree" ]] && "$tree"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --aln "$infile" \
    $tree_arg \
    $clock_rate_arg \
    --dates "$date_file" \
    --outdir "$outdir" \
    --confidence \
    --covariation \
    $more_args 

# Report output
log_time "Showing the contents of the output dates.tsv:"
cat -n "$outdir"/dates.tsv

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
