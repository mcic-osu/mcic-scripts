#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=align
#SBATCH --output=slurm-align-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Align nucleotide or amino acid sequences in a multi-FASTA file with MAFFT (default) or MUSCLE"
SCRIPT_VERSION="2025-02-16"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
container_dir="$HOME/containers"
version_only=false                 # When true, just print tool & script version info and exit

# Defaults - tool parameters
aligner=mafft
fix_header=true

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
    echo "      sbatch $0 -i my.fa -o aln.fa"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input multi-FASTA file with sequences to be aligned"
    echo "                              FASTA can contain either nucleotide or amino acid (protein) sequences"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)  [default: same as input dir]"
    echo "  --aligner           <str>   Aligner: 'mafft' or 'muscle'            [default: $aligner]"
    echo "  --no_header_fix             Don't fix output FASTA header           [default: fix header]"
    echo "                              By default, the script will remove"
    echo "                              aligner-added header text like 'reverse'"
    echo "  --more_opts         <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
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
more_opts=
threads=
container_path=
container_url=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --aligner )         shift && aligner=$1 ;;
        --no_header_fix )   shift && fix_header=false ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env )             shift && env=$1 ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 ;;
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

# Constants and vars based on aligner choice
TOOL_BINARY="$aligner"
if [[ "$aligner" == "mafft" ]]; then
    TOOL_NAME=MAFFT
    VERSION_COMMAND="$TOOL_BINARY --version"
    conda_path=/fs/ess/PAS0471/jelmer/conda/mafft
elif [[ "$aligner" == "muscle" ]]; then
    TOOL_NAME=MUSCLE
    VERSION_COMMAND="$TOOL_BINARY --version"
    conda_path=/fs/ess/PAS0471/jelmer/conda/muscle
else
    die "Aligner should be 'mafft' or 'muscle' but is $aligner"
fi

# Load software
load_env "$conda_path" "$container_path" #"$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# Define outputs based on script parameters
[[ -z "$outdir" ]] && outdir=$(dirname "$(realpath "$infile")")
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
outfile="$outdir"/$(basename "${infile%.*}")_aln.fa

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input file:                               $infile"
echo "Output file:                              $outfile"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$infile"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the alignment
log_time "Running $TOOL_NAME..."
if [[ "$aligner" == "mafft" ]]; then
    runstats $CONTAINER_PREFIX $TOOL_BINARY \
        --reorder \
        --auto \
        --adjustdirection \
        --leavegappyregion \
        --thread "$threads" \
        $more_opts \
        "$infile" > "$outfile"
elif [[ "$aligner" == "muscle" ]]; then
    runstats $CONTAINER_PREFIX $TOOL_BINARY \
        -align "$infile" \
        -output "$outfile"
        $more_opts
fi

# Remove aligner-added extra info after space from FASTA header lines
# ... and remove "_R_" prefixes for reverse-complemented seqs
if [[ "$fix_header" == true ]]; then
    log_time "Fixing FASTA headers..."
    sed -i -E -e 's/(^>[^ ]+) .*/\1/' -e 's/^>_R_/>/' "$outfile"
fi

# Final reporting
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
