#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=pycoqc
#SBATCH --output=slurm-pycoqc-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run PycoQC to QC an ONT sequencing summary file"
SCRIPT_VERSION="2024-06-29"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=pycoQC
TOOL_NAME=pycoQC
TOOL_DOCS=https://a-slide.github.io/pycoQC
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env_type=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/pycoqc
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                 # When true, just print tool & script version info and exit

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
    echo "      sbatch $0 -i data/sequencing_summary.txt -o results/pycoqc"
    echo "  - Only reads with a Phred quality score of 10 or more will be marked as PASS:"
    echo "      sbatch $0 -i data/sequencing_summary.txt -o results/pycoqc --min_qual 10"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input file: a 'sequencing summary' text file output by the sequencing machine"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --min_qual          <int>   Min. Phred read quality score to PASS   [default: not used]"
    echo "  --min_len           <int>   Min. read length to PASS                [default: not used]"
    echo "  --more_opts         <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env_type               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env_type]"
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
infile=
outdir=
min_len= && len_arg=
min_qual= && qual_arg=
more_opts=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --min_qual )        shift && min_qual=$1 ;;
        --min_len )         shift && min_len=$1 ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env_type )             shift && env_type=$1 ;;
        --no_strict )       strict_bash=false ;;
        --dl_container )    dl_container=true ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 && dl_container=true ;;
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
[[ "$strict_bash" == true ]] && set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# Determine the output file name
infile_noext=${infile%.*}
outfile="$outdir"/$(basename "$infile_noext")_pycoqc.html

# Other options
[[ -n "$min_len" ]] && len_arg="--min_pass_len $min_len"
[[ -n "$min_qual" ]] && qual_arg="--min_pass_qual $min_qual"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input file:                               $infile"
echo "Output dir:                               $outdir"
[[ -n $min_len ]] && echo "Min. read length to pass:                 $min_len"
[[ -n $min_qual ]] && echo "Min. read quality to pass:                $min_qual"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    -f "$infile" \
    -o "$outfile" \
    $qual_arg \
    $len_arg \
    $more_opts

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
