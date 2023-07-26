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
SCRIPT_VERSION="2023-07-25"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=treetime
TOOL_NAME=Treetime
TOOL_DOCS=https://treetime.readthedocs.io
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=container                          # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/project/PAS0471/jelmer/conda/treetime
container_path=/fs/ess/PAS0471/containers/treetime_0.10.1--pyh7cba7a3_0.sif
container_url=docker://quay.io/biocontainers/treetime:0.10.1--pyh7cba7a3_0
dl_container=false
container_dir="$HOME/containers"

# Defaults - tool parameters
#TODO

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo "                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage:"
    echo "      sbatch $0 -i results/snippy/core.full.aln -o results/treetime --dates metadata/dates.tsv"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--alignment      <file>  Input alignment file in FASTA, Phylip, or VCF format"
    echo "  --dates             <file>  Input CSV/TSV with sampling date for each sample"
    echo "                                Should have columns 'node_name' and 'date', with dates as YYYY-MM-DD"
    echo "                                If column names are different, use '--opts --date-column <date-column> --name-column <samplename-column>'"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --tree              <file>  Tree file in Nexus or Newick format"
    echo "                                If no tree file is provided, TreeTime will infer one"
    echo "  --clock_rate        <str>   Clock rate to assume                        [default: no clock rate]"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or --dl_container is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: false]"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
    echo
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
outdir=
date_file=
tree= && tree_arg=
clock_rate= && clock_rate_arg=
opts=
version_only=false

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --dates )           shift && date_file=$1 ;;
        --tree )            shift && tree=$1 ;;
        --clock_rate )      shift && clock_rate=$1 ;;
        --opts )            shift && opts=$1 ;;
        --env )             shift && env=$1 ;;
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
set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input alignment file specified, do so with -i/--infile" "$all_args"
[[ -z "$date_file" ]] && die "No input dates file specified, do so with --dates" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -f "$date_file" ]] && die "Input file $date_file does not exist"
[[ -n "$tree" && ! -f "$tree" ]] && die "Input tree file $tree does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
[[ -n $tree ]] && tree_arg="--tree $tree"
[[ -n $clock_rate ]] && clock_rate_arg="--clock-rate $clock_rate"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input alignment file:                     $infile"
echo "CSV/TSV dates input file:                 $date_file"
[[ -n "$tree" ]] && echo "Input tree file:                          $tree"
echo "Output dir:                               $outdir"
[[ -n "$clock_rate" ]] && echo "Clock rate:                               $clock_rate"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:        $opts"
log_time "Listing the input file(s):"
ls -lh "$infile" "$date_file"
[[ -n "$tree" ]] && ls -lh "$tree"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    --aln "$infile" \
    $tree_arg \
    $clock_rate_arg \
    --dates "$date_file" \
    --outdir "$outdir" \
    --confidence \
    --covariation \
    $opts

#? TreeTime options:
#> --reroot         min_dev / least-squares / oldest
#> --keep-root

log_time "Showing the contents of the output dates.tsv:"
cat -n "$outdir"/dates.tsv

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
