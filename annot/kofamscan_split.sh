#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=kofamscan_split
#SBATCH --output=slurm-kofamscan_split-%j.out

# Run KoFamScan on a protein FASTA to assign KEGG K-numbers to the proteins

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=kofamscan_split.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/seqkit
readonly TOOL_BINARY=seqkit
readonly TOOL_NAME=KoFamScan

readonly KOFAMSCAN_SCRIPT=mcic-scripts/annot/kofamscan.sh

# Option defaults
db_dir=/fs/ess/PAS0471/jelmer/refdata/kegg
download_db=false
out_format=mapper

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  Run KoFamScan on a protein FASTA file to assign KEGG K-numbers to the proteins"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input-file> -o <output-dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input protein FASTA file"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --db_dir        <dir>   Dir with (or for) the KEGG database         [default: /fs/ess/PAS0471/jelmer/refdata/kegg]"
    echo "  --download_db           Download the KEGG database"
    echo "  --out_format    <str>   Output format: 'detail', 'detail-tsv', or 'mapper' [default: 'mapper']"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly/proteins.faa -o results/kofamscan/ko.tsv"
    echo
    echo "OUTPUT:"
    echo "  - ..."
}

# Load software
load_tool_conda() {
    set +u
    module load miniconda3/4.12.0-py39 # Load the OSC Conda module
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
    echo "Run using $SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION (https://github.com/mcic-osu/mcic-scripts)"
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
outdir=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && readonly infile=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        --batch_size )      shift && readonly batch_size=$1 ;;
        --db_dir )          shift && readonly db_dir=$1 ;;
        --download_db )     readonly download_db=true ;;
        --out_format )      shift && readonly out_format=$1 ;;
        -v )                script_version; exit 0 ;;
        -h )                script_help; exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ -z "$download_db" && ! -d "$db_dir" ]] && die "Database dir $db_dir does not exist, use --download_db if you need to download it"

# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
load_tool_conda

# ==============================================================================
#                      DEFINE OUTPUTS AND DERIVED INPUTS
# ==============================================================================
# Derived inputs
readonly profile_dir="$db_dir"/profiles
readonly ko_list="$db_dir"/ko_list
[[ ! -d "$profile_dir" ]] && die "Input dir $profile_dir does not exist"
[[ ! -f "$ko_list" ]] && die "Input file $ko_list does not exist"

# Get the file extension
file_ext=${infile##*.}

# Define outputs based on script parameters
readonly log_dir="$outdir"/logs

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input file:                               $infile"
echo "Output dir:                               $outdir"
echo "Output format:                            $out_format"
echo "Database dir:                             $db_dir"
echo "Download the database?                    $download_db"
echo "FASTA batch size:                         $batch_size"
echo
echo "# Listing the input file(s):"
ls -lh "$infile"

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir" "$outdir"/split_fasta

# Download the databases
if [[ "$download_db" == true ]]; then
    log_time "Downloading the database..."
    mkdir -pv "$db_dir"
    cd "$db_dir" || die "Can't change to $db_dir"
    
    wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    gunzip ko_list.gz
    
    wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    tar xf profiles.tar.gz
    
    cd -
fi

# Split the FASTA file
log_time "Splitting the FASTA files into chunks of $batch_size sequences..."
runstats "$TOOL_BINARY" split2 \
    "$infile" \
    -s "$batch_size" \
    --out-dir "$outdir"/split_fasta

log_time "Created $(ls "$outdir"/split_fasta/*."$file_ext" | wc -l) FASTA files"

# Run the tool
for fa in "$outdir"/split_fasta/*."$file_ext"; do
    file_id=$(basename "$fa" ."$file_ext")
    outfile="$outdir"/"$file_id"/"$file_id".tsv
    
    log_time "Submitting script for file $fa..."
    sbatch "$KOFAMSCAN_SCRIPT" \
        -i "$fa" \
        -o "$outfile" \
        --db_dir "$db_dir" \
        --out_format "$out_format"
done


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Done with script $SCRIPT_NAME"
echo
