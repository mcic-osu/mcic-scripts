#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=170G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=racon
#SBATCH --output=slurm-racon-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Racon (Minimap then 1 or more rounds of Racon) to polish a genome
assembly either with short or long reads"
SCRIPT_VERSION="2023-09-26"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=racon
TOOL_NAME=Racon
TOOL_DOCS=https://github.com/lbcb-sci/racon
VERSION_COMMAND="$TOOL_BINARY --version; echo '# Version of Minimap:'; minimap2 --version"

# Defaults - generics
env_type=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/racon
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
version_only=false                 # When true, just print tool & script version info and exit

# Defaults - tool parameters
minimap_preset="map-ont"
iterations=2

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
    echo "      batch $0 --assembly results/flye/assembly.fasta --reads data/fastq/my.fastq.gz -o results/racon"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --assembly          <file>  Input assembly: FASTA file (to be corrected)"
    echo "  --reads             <file>  Input reads: FASTQ file (reads used for correction)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --iterations        <int>   Number of Racon iterations (1 or 2)     [default: 2]"
    echo "  --minimap_preset    <str>   Minimap preset                          [default: 'map-ont']"
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
reads=
assembly=
outdir=
more_opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --reads )           shift && reads=$1 ;;
        --assembly )        shift && assembly_in=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --minimap_preset )  shift && minimap_preset=$1 ;;
        --iterations )      shift && iterations=$1 ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env_type )             shift && env_type=$1 ;;
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
set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$assembly_in" ]] && die "No input assembly FASTA file specified, do so with --assembly" "$all_opts"
[[ -z "$reads" ]] && die "No input FASTQ file specified, do so with --reads" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$assembly_in" ]] && die "Input assembly FASTA file $assembly_in does not exist"
[[ ! -f "$reads" ]] && die "Input reads FASTQ file $reads does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR" "$outdir"/minimap
assembly_ext=$(echo "$assembly_in" | sed -E 's/.*(\.fn?a?s?t?a$)/\1/')
assembly_id=$(basename "$assembly_in" "$assembly_ext")
assembly_out1="$outdir"/"$assembly_id"_racon1.fasta
[[ "$iterations" -eq 2 ]] && assembly_out2="$outdir"/"$assembly_id"_racon2.fasta
align_1="$outdir"/minimap/"$assembly_id"_iter1.sam
[[ "$iterations" -eq 2 ]] && align_2="$outdir"/minimap/"$assembly_id"_iter2.sam
[[ "$iterations" -gt 2 ]] && die "Number of Racon iterations cannot be greater than 2 (You asked for $iterations)" 

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input reads (FASTQ) file:                 $reads"
echo "Input assembly (FASTA) file:              $assembly_in"
echo "Output dir:                               $outdir"
echo "Nr of Racon iterations:                   $iterations"
echo "Minimap preset:                           $minimap_preset"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$reads" "$assembly_in"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               FUNCTIONS
# ==============================================================================
# Function to run Racon
Run_racon() {
    assembly_in=${1:-none}
    alignments=${2:-none}
    assembly_out=${3:-none}

    [[ $assembly_in == "none" ]] && die "No assembly for function Run_racon"
    [[ $alignments == "none" ]] && die "No alignments for function Run_racon"
    [[ $assembly_out == "none" ]] && die "No outfile for function Run_racon"

    runstats $TOOL_BINARY \
        "$reads" \
        "$alignments" \
        "$assembly_in" \
        -t "$threads" \
        $more_opts \
        > "$assembly_out"
}

# Function un Minimap
Run_minimap() {
    assembly=${1:-none}
    align_out=${2:-none}

    [[ $assembly == "none" ]] && die "No assembly for function Run_minimap"
    [[ $align_out == "none" ]] && die "No outfile for function Run_minimap"

    runstats minimap2 \
        -x "$minimap_preset" \
        -t "$threads" \
        -a \
        "$assembly" \
        "$reads" \
        > "$align_out"
}

# ==============================================================================
#                               RUN
# ==============================================================================
# Minimap iteration 1
if [[ ! -s "$align_1" ]]; then
    log_time "Now running the first iteration of Minimap..."
    Run_minimap "$assembly_in" "$align_1"
else
    log_time "Minimap SAM from iteration 1 exists, skipping step..."
    ls -lh "$align_1"
fi

# Racon iteration 1
if [[ ! -s "$assembly_out1" ]]; then
    echo -e "\n====================================================================="
    log_time "Now running the first iteration of Racon..."
    Run_racon "$assembly_in" "$align_1" "$assembly_out1"
else
    log_time "Assembly from Racon iteration 1 exists, skipping step..."
    ls -lh "$assembly_out1"
fi

if [[ "$iterations" -eq 2 ]]; then
    # Minimap iteration 2
    if [[ ! -s "$align_2" ]]; then
        echo -e "\n====================================================================="
        log_time "Now running the second iteration of Minimap..."
        Run_minimap "$assembly_out1" "$align_2"
    else
        log_time "Minimap SAM from iteration 2 exists, skipping step..."
        ls -lh "$align_2"
    fi

    # Racon iteration 2
    if [[ ! -s "$assembly_out2" ]]; then
        echo -e "\n====================================================================="
        log_time "Now running the second iteration of Racon..."
        Run_racon "$assembly_out1" "$align_2" "$assembly_out2"
    else
        log_time "Assembly from Racon iteration 2 exists, skipping step..."
        ls -lh "$assembly_out2"
    fi
fi

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
