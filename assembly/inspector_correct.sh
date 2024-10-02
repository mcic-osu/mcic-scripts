#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=inspector_corr
#SBATCH --output=slurm-inspector_corr-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Inspector-correct to correct a genome assembly"
SCRIPT_VERSION="2024-09-29"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=inspector-correct.py
TOOL_NAME=Inspector
TOOL_DOCS=https://github.com/Maggi-Chen/Inspector
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
#? Using a container, got a memory error with the Conda env
env=container                       # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/inspector
container_path=
container_url=oras://community.wave.seqera.io/library/inspector:1.3.1--1e8dfa3d2ec456ac
dl_container=false
container_dir="$HOME/containers"
version_only=false                 # When true, just print tool & script version info and exit

# Defaults for tool
base_error=false

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
    echo "      sbatch $0 -i results/inspector -o results/inspector/asm.fasta --datatype 'nano-corr'"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--inspector_dir  <dir>   Dir with Inspector results -- run inspector.sh before this script"
    echo "  -o/--outfile        <dir>   Output assembly FASTA file."
    echo "                              NOTE: Just provide the filename, no path: file will be added to Inspector dir"
    echo "  --datatype          <str>   Input read type: pacbio-raw, pacbio-hifi, pacbio-corr, nano-raw, nano-corr [default: 'nano-raw']"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --base_error                Also correct base-errors                [default: don't correct]"
    echo "  --more_opts         <str>   Quoted string with additional options for $TOOL_NAME"
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
indir=
outfile=
base_error_opt=
more_opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --inspector_dir )  shift && indir=$1 ;;
        -o | --outfile )        shift && outfile=$1 ;;
        --datatype )            shift && datatype=$1 ;;
        --base_error )          base_error=true ;;
        --more_opts )           shift && more_opts=$1 ;;
        --env )                 shift && env=$1 ;;
        --dl_container )        dl_container=true ;;
        --container_dir )       shift && container_dir=$1 ;;
        --container_url )       shift && container_url=$1 && dl_container=true ;;
        -h | --help )           script_help; exit 0 ;;
        -v )                    script_version; exit 0 ;;
        --version )             version_only=true ;;
        * )                     die "Invalid option $1" "$all_opts" ;;
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
[[ -z "$indir" ]] && die "No input dir specified, do so with -i/--inspector_dir" "$all_opts"
[[ -z "$outfile" ]] && die "No output file specified, do so with -o/--outfile" "$all_opts"
[[ ! -d "$indir" ]] && die "Input dir $indir does not exist"

# Define outputs based on script parameters
#outdir=$(dirname "$outfile")
indir=$(realpath "$indir")
outdir="$indir"
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
[[ "$base_error" == false ]] && base_error_opt="--skip_baseerror"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input Inspector results dir:              $indir"
echo "Output assembly file:                     $outfile"
echo "Data type:                                $datatype"
echo "Correct base-errors, too:                 $base_error"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$indir"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Move to the outdir or some files will go to the working dir
log_time "Moving into the outdir..."
cd "$outdir" || exit 1

log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    --inspector . \
    -o . \
    --datatype "$datatype" \
    --thread "$threads" \
    "$base_error_opt" \
    $more_opts

log_time "Renaming the output file..."
mv -v contig_corrected.fa "$outfile"

log_time "Listing files in the output dir:"
ls -lhd $PWD/*
final_reporting "$LOG_DIR"
