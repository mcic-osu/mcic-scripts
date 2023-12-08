#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=dl_refseq
#SBATCH --output=slurm-dl_refseq-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Download and concatenate the NCBI RefSeq database"
SCRIPT_VERSION="2023-12-06"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh

# Defaults - generics
strict_bash=true

# Constants - tool parameters
BASE_URL=ftp://ftp.ncbi.nlm.nih.gov/refseq/release

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
    echo "      sbatch $0 --lib xx -o results/refseq"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --lib               <str>   RefSeq library"
    echo "                                Options are the dir names here: https://ftp.ncbi.nlm.nih.gov/refseq/release/"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --no_strict                 Don't use strict Bash settings ('set -euo pipefail') -- can be useful for troubleshooting"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
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
refseq_lib=
outdir=
more_opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --lib )             shift && refseq_lib=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --more_opts )       shift && more_opts=$1 ;;
        --no_strict )       strict_bash=false ;;
        -h | --help )       script_help; exit 0 ;;
        -v )                script_version; exit 0 ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
[[ "$strict_bash" == true ]] && set -euo pipefail

# Check options provided to the script
[[ -z "$refseq_lib" ]] && die "No RefSeq library specified, do so with --lib" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
outfile="$outdir"/refseq_"$refseq_lib".fa

# NCBI key
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
#? https://ftp.ncbi.nlm.nih.gov/refseq/release/
#? See also https://github.com/DerrickWood/kraken2/issues/115

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "RefSeq library:                           $refseq_lib"
echo "Output dir:                               $outdir"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
log_time "Starting the download..."
runstats wget \
    "$BASE_URL"/"$refseq_lib"/*genomic.fna.gz \
    -P "$outdir"

# Concatenate the FASTQ files
log_time "Concatenating the FASTA files..."
runstats cat "$outdir"/*genomic.fna.gz > "$outfile"

log_time "Listing the main output file:"
ls -lh "$outfile"
final_reporting "$LOG_DIR"
