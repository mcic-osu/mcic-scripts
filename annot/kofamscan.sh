#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=kofamscan
#SBATCH --output=slurm-kofamscan-%j.out

# Run KoFamScan on a protein FASTA to assign KEGG K-numbers to the proteins

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=kofamscan.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/kofamscan
readonly TOOL_BINARY=exec_annotation
readonly TOOL_NAME=KoFamScan
readonly TOOL_DOCS=https://github.com/takaram/kofam_scan
readonly TOOL_PAPER=https://academic.oup.com/bioinformatics/article/36/7/2251/5631907

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
    echo "  -o/--outfile    <dir>   Output file"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --db_dir        <dir>   Dir with (or for) the KEGG database         [default: /fs/ess/PAS0471/jelmer/refdata/kegg]"
    echo "  --download_db           Download the KEGG database"
    echo "  --out_format    <str>   Output format: 'detail', 'detail-tsv', or 'mapper' [default: 'mapper']"
    echo "  --more_args     <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for $TOOL_NAME and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/assembly/proteins.faa -o results/kofamscan/ko.tsv"
    echo
    echo "OUTPUT:"
    echo "  - A single tabular output file with gene-to-KO mappings,"
    echo "    whose format depends on the --out_format option."
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
    echo "  - See also: https://taylorreiter.github.io/2019-05-11-kofamscan/"
    echo
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

# Print the tool's version
print_version() {
    set +e
    load_tool_conda
    "$TOOL_BINARY" -v
    set -e
}

# Print the tool's help
tool_help() {
    load_tool_conda
    "$TOOL_BINARY" -h
}

# Print SLURM job resource usage info
resource_usage() {
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
}

# Print SLURM job requested resources
slurm_resources() {
    set +u
    log_time "SLURM job information:"
    echo "Account (project):                        $SLURM_JOB_ACCOUNT"
    echo "Job ID:                                   $SLURM_JOB_ID"
    echo "Job name:                                 $SLURM_JOB_NAME"
    echo "Memory (MB per node):                     $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):                          $SLURM_CPUS_PER_TASK"
    echo "Time limit:                               $SLURM_TIMELIMIT"
    echo -e "=================================================================\n"
    set -u
}

# Set the number of threads/CPUs
set_threads() {
    set +u
    if [[ "$is_slurm" == true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            readonly threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            readonly threads="$SLURM_NTASKS"
        else 
            log_time "WARNING: Can't detect nr of threads, setting to 1"
            readonly threads=1
        fi
    else
        readonly threads=1
    fi
    set -u
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
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outfile )    shift && outfile=$1 ;;
        --db_dir )          shift && db_dir=$1 ;;
        --download_db )     readonly download_db=true ;;
        --out_format )      shift && readonly out_format=$1 ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -h )                script_help; exit 0 ;;
        -v | --version )         print_version; exit 0 ;;
        --help )            tool_help; exit 0;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$outfile" ]] && die "No output file specified, do so with -o/--outfile" "$all_args"
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
set_threads

# ==============================================================================
#                      DEFINE OUTPUTS AND DERIVED INPUTS
# ==============================================================================
# Derived inputs
profile_dir="$db_dir"/profiles
ko_list="$db_dir"/ko_list
[[ ! -d "$profile_dir" ]] && die "Input dir $profile_dir does not exist"
[[ ! -f "$ko_list" ]] && die "Input file $ko_list does not exist"

# Make paths absolute
infile=$(realpath "$infile")
[[ ! "$outfile" =~ ^/ ]] && outfile="$PWD"/"$outfile"
[[ ! "$db_dir" =~ ^/ ]] && db_dir="$PWD"/"$db_dir"

# Define outputs based on script parameters
outdir=$(dirname "$outfile")
readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input file:                               $infile"
echo "Nr of entries in the input file:          $(grep -c "^>" "$infile")"
echo "Output file:                              $outfile"
echo "Output format:                            $out_format"
echo "Database dir:                             $db_dir"
echo "Download the database?                    $download_db"
[[ $more_args != "" ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
echo "Number of threads/cores:                  $threads"
echo
echo "# Listing the input file(s):"
ls -lh "$infile"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir"

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

# Copy the DB to the TMPDIR
if [[ "$is_slurm" == true ]]; then
    log_time "Copying the database to the TMPDIR $TMPDIR..."
    cp -r "$db_dir" "$TMPDIR"
    db_dir="$TMPDIR"/$(basename "$db_dir")

    profile_dir="$db_dir"/profiles
    ko_list="$db_dir"/ko_list
    [[ ! -d "$profile_dir" ]] && die "Input dir $profile_dir does not exist"
    [[ ! -f "$ko_list" ]] && die "Input file $ko_list does not exist"
fi

# Move into the outdir
# This is needed because kofamscan will create a 'tmp' dir in the working dir
if [[ "$is_slurm" == true ]]; then
    cd "$TMPDIR" || Die "Can't move to $TMPDIR"
else
    cd "$outdir" || Die "Can't move to $outdir"
fi

# Run the tool
log_time "Running $TOOL_NAME..."
runstats "$TOOL_BINARY" \
    -o "$outfile" \
    --profile="$profile_dir" \
    --ko-list="$ko_list" \
    --cpu="$threads" \
    --format="$out_format" \
    $more_args \
    "$infile"

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
print_version | tee "$version_file"
script_version | tee -a "$version_file" 
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME"
echo
