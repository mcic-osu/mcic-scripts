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
SCRIPT_VERSION="2025-03-05"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh

# Defaults - generics
env_type=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/mafft # Also contains Muscle
container_dir="$HOME/containers"
container_url=

# Defaults - tool parameters
aligner=mafft
fix_header=true

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "
                        $0
    v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)
            =================================================

DESCRIPTION:
$DESCRIPTION
    
USAGE / EXAMPLE COMMANDS:
  - Basic usage example:
      sbatch $0 -i TODO -o results/TODO
    
REQUIRED OPTIONS:
  -i/--infile         <file>  Input multi-FASTA file with sequences to be aligned
                              FASTA can contain either nucleotide or amino acid (protein) sequences
    
OTHER KEY OPTIONS:
  -o/--outdir         <dir>   Output dir (will be created if needed)            [default: same as input dir]
  --aligner           <str>   Aligner: 'mafft' or 'muscle'                      [default: $aligner]
  --no_header_fix             Don't fix output FASTA header                     [default: fix header]
                              By default, the script will remove aligner-added
                              header text like 'reverse'

UTILITY OPTIONS:
  --env_type          <str>   Use a Singularity container ('container')         [default: $env_type]
                              or a Conda environment ('conda') 
  --conda_path        <dir>   Full path to a Conda environment to use           [default: $conda_path]
  --container_url     <str>   URL to download a container from                  [default: $container_url]
  --container_dir     <str>   Dir to download the container to                  [default: $container_dir]
  --container_path    <file>  Singularity container image file (.sif) to use
  -h/--help                   Print help and exit
  -v/--version                Print the version of this script and of $TOOL_NAME
    
TOOL DOCUMENTATION:
$TOOL_DOCS
"
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
version_only=false   # When true, just print tool & script version info and exit
infile=
outdir=
more_opts=
threads=
container_path=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --aligner )         shift && aligner=$1 ;;
        --no_header_fix )   shift && fix_header=false ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env_type )        shift && env_type=$1 ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 ;;
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

# Constants and vars based on aligner choice
TOOL_BINARY="$aligner"
if [[ "$aligner" == "mafft" ]]; then
    TOOL_NAME=MAFFT
elif [[ "$aligner" == "muscle" ]]; then
    TOOL_NAME=MUSCLE
else
    die "Aligner should be 'mafft' or 'muscle' but is $aligner"
fi
VERSION_COMMAND="$TOOL_BINARY --version"

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
