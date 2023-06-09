#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=longest_transcripts.sh
#SBATCH --output=slurm-longest_transcripts.sh-%j.out

# From a proteome (or nucleotide) FASTA file,
# create a new FASTA file with only the longest transcript (isoform) per gene

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=longest_transcripts.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/seqkit
readonly TOOL_BINARY=seqkit

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  From a proteome (or nucleotide) FASTA file,"
    echo "  create a new FASTA file with only the longest transcript (isoform) per gene"
    echo "  This will work on FASTA files produced by Braker, with sequence IDs in"
    echo "  the format of 'g1.t1' (gene 1, transcript 1) etc"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i <input-file> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input FASTA file"
    echo "                          The sequence IDs are assumed to be in the format of 'g1.t1' (gene 1, transcript 1) etc"
    echo "  -o/--outfile    <dir>   Output FASTA file (dir will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "OUTPUT:"
    echo "  - A FASTA file with only the longest transcript for each gene"
    echo "  - A TSV file with the lengths for each of the sequences in the input file" 
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
    echo
}

# Load software
load_tool_conda() {
    set +u
    module load "$MODULE" # Load the OSC Conda module
    # Deactivate any active Conda environments:
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi
    source activate "$CONDA_ENV" # Activate the focal environment
    set -u
}

# Exit upon error with a message
die() {
    local error_message=${1}
    local error_args=${2-none}
    log_time "$0: ERROR: $error_message" >&2
    log_time "For help, run this script with the '-h' option" >&2
    if [[ "$error_args" != "none" ]]; then
        log_time "Arguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    log_time "EXITING..." >&2
    exit 1
}

# Log messages that include the time
log_time() { echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')]" ${1-""}; }

# Print the script version
script_version() {
    echo "Run using $SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION ($SCRIPT_URL)"
}

# Print the tool's version
tool_version() {
    set +e
    load_tool_conda
    "$TOOL_BINARY" version
    set -e
}

# Print SLURM job resource usage info
resource_usage() {
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
}

# Resource usage information for any process
runstats() {
    /usr/bin/time -f \
        "\n# Ran the command: \n%C
        \n# Run stats by /usr/bin/time:
        Time: %E   CPU: %P    Max mem: %M K    Exit status: %x \n" \
        "$@"
}

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
infile=
outfile=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && readonly infile=$1 ;;
        -o | --outfile )    shift && readonly outfile=$1 ;;
        -v )                script_version; exit 0 ;;
        -h )                script_help; exit 0 ;;
        --version )         tool_version; exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$outfile" ]] && die "No output dir specified, do so with -o/--outfile" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Load software
load_tool_conda

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Define outputs based on script parameters
outdir=$(dirname "$outfile")
seqlen_file="$outdir"/seqlengths.tsv
longtrans_file="$outdir"/longest_transcripts.txt

readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input FASTA file:                         $infile"
echo "Output FASTA file:                        $outfile"
echo "Output TSV file with sequence lengths:    $seqlen_file"
echo "Output text file with longest transcript IDS: $longtrans_file"
echo
log_time "Listing the input file(s):"
ls -lh "$infile"

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir"

# Run the tool
log_time "Getting the sequence lengths..."
runstats $TOOL_BINARY fx2tab --length --name "$infile" > "$seqlen_file"

log_time "Getting the IDs of the longest transcripts..."
tr '.' '\t' < "$seqlen_file" |
    sort -k1,1 -k 3,3nr |
    sort -k1,1 -u |
    cut -f1,2 |
    tr '\t' '.' \
    > "$longtrans_file"

log_time "Extracting the longest transcripts..."
$TOOL_BINARY grep -f "$longtrans_file" "$infile" > "$outfile"

# Report and check
n_longest=$(wc -l < "$longtrans_file")
n_out=$(grep -c ">" "$outfile")
[[ "$n_longest" != "$n_out" ]] && die "Output FASTA does not contain the same nr of sequences ($n_out) as the file with longest transcripts ($n_longest)"

echo
echo "Number of sequences in the input file:      $(grep -c ">" "$infile")"
echo "Number of sequences in the output file:     $(grep -c ">" "$outfile")"

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
tool_version | tee "$version_file"
script_version | tee -a "$version_file" 
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME"
echo
