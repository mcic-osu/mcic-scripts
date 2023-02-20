#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=resfinder
#SBATCH --output=slurm-resfinder-%j.out

# Run ResFinder to detect AMR genes and point mutations in a bacterial genome assembly

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME="resfinder.sh"
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly CONDA_ENV="/fs/ess/PAS0471/jelmer/conda/resfinder-4.1.11"
readonly TOOL_BINARY="run_resfinder.py"
readonly TOOL_NAME="ResFinder"
readonly TOOL_DOCS="https://bitbucket.org/genomicepidemiology/resfinder"
readonly TOOL_PAPER=""

# Option defaults
get_db=true            # Don't download the database
point_mut=true         # Include point mutations
acquired=true
min_cov=0.90           # Coverage threshold for BLAST hits
min_id=0.90            # Identity threshold for BLAST hits

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                      Run $TOOL_NAME"
    echo "======================================================================"
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input file: a nucleotide FASTA file with a genome assembly"
    echo "  -o/--outdir     <dir>   Output dir (use a separate dir per assembly)"
    echo "  --species       <str>   Species name -- see ResFinder docs for possible values"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --no_point_mut          Skip point-mutation-based resistance    [default: include]"
    echo "  --no_aqcuired           Skip aqcuired genes                     [default: include]"
    echo "  --min_cov       <num>   Coverage threshold                      [default: 0.90]"
    echo "  --min_id        <num>   Identity threshold                      [default: 0.90]"
    echo "  --db_dir                Dir with/for ResFinder DBs              [default: <outdir>/dbs]"
    echo "  --dont_get_db           Don't download the ResFinder DBs        [default: download]"
    echo "  --more_args     <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for $TOOL_NAME and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  $0 -i results/spades/assembly.fa -o results/resfinder"
    echo
    echo "  for asm in result/assemblies/*fasta; do"
    echo "      outdir=results/resfinder/$(basename "$asm" .fasta)"
    echo "      sbatch mcic-scripts/bact/resfinder.sh -i "$asm" -o "$outdir""
    echo "  done"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
    echo
}

# Load software
load_tool_conda() {
    set +u
    # Load the OSC Conda module 
    module load miniconda3/4.12.0-py39
    # Deactivate any activae Conda environments
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi
    # Activate the focal environment
    source activate "$CONDA_ENV"
    set -u
}

# Exit upon error with a message
die() {
    local error_message=${1}
    local error_args=${2-none}
    
    echo -e "\n============================================================" >&2
    log_time "$0: ERROR: $error_message" >&2
    log_time "For help, run this script with the '-h' option" >&2
    if [[ "$error_args" != "none" ]]; then
        log_time "Arguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    log_time "EXITING..." >&2
    echo -e "============================================================\n" >&2
    exit 1
}

# Log messages that include the time
log_time() {
    echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')]" ${1-""}
}

# Print the script version
script_version() {
    echo "Run using $SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION (https://github.com/mcic-osu/mcic-scripts)"
}

# Print the tool's version
tool_version() {
    set +e
    load_tool_conda
    conda list | grep -i $TOOL_NAME | tail -n +2
    set -e
}

# Print the tool's help
tool_help() {
    load_tool_conda
    "$TOOL_BINARY" --help
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
outdir=
species=
db_dir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --species )         shift && species=$1 ;;
        --dont_get_db )     get_db=false ;;
        --db_dir )          shift && db_dir=$1 ;;
        --no_point_mut )    point_mut=false && point_arg="" ;;
        --no_aqcuired )     acquired=false && acquired_arg="" ;;
        --min_cov )         shift && readonly min_cov=$1 ;;
        --min_id )          shift && readonly min_id=$1 ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -v )                script_version; exit 0 ;;
        -h )                script_help; exit 0 ;;
        --version )         tool_version; exit 0 ;;
        --help )            tool_help; exit 0;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ -z "$species" ]] && die "No species specified, do so with --species" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ -n "$db_dir" && ! -d "$db_dir" ]] && die "DB dir $db_dir does not exist"
[[ "$get_db" == false && -z "$db_dir" ]] && die "You have to specify a DB dir with --db_dir when using --dont_get_db"

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
#                   DEFINE OUTPUTS, DERIVED INPUTS, AND ARGS
# ==============================================================================
# Define outputs based on script parameters
[[ -z "$db_dir" ]] && db_dir="$outdir"/dbs

# Make paths absolute
infile=$(realpath "$infile")
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"
[[ ! "$db_dir" =~ ^/ ]] && db_dir="$PWD"/"$db_dir"

# Define outputs based on script parameters
readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs

# Define args
[[ "$point_mut" == true ]] && point_arg="--point"
[[ "$acquired" == true ]] && acquired_arg="--acquired"

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input genome assembly FASTA file:         $infile"
echo "Output dir:                               $outdir"
echo "Species:                                  $species"
echo "Include point mutations:                  $point_mut"
echo "Include acquired genes:                   $acquired"
echo "Min. coverage threshold:                  $min_cov"
echo "Min. identity threshold:                  $min_id"
[[ $more_args != "" ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
echo "Number of threads/cores:                  $threads"
echo
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir" "$db_dir"

# Download/update the database
if [[ "$get_db" == true ]]; then
    log_time "Downloading the databases..."

    cd "$db_dir" || die "Can't move to DB dir $db_dir"
    
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db
    cd resfinder_db || exit 1
    python3 INSTALL.py

    if [[ "$point_mut" == true ]]; then
        cd "$db_dir" || die "Can't move to DB dir $db_dir"
        git clone https://bitbucket.org/genomicepidemiology/pointfinder_db
        cd pointfinder_db || exit 1
        python3 INSTALL.py
    fi
fi

# Move the output dir
cd "$outdir" || die "Can't move to output dir $outdir"

# Run
log_time "Running $TOOL_NAME..."
runstats run_resfinder.py \
    --inputfasta "$infile" \
    --outputPath "$outdir" \
    --species "$species" \
    --min_cov "$min_cov" \
    --threshold "$min_id" \
    --db_path_res "$db_dir"/resfinder_db \
    --db_path_point "$db_dir"/dbs/pointfinder_db \
    $acquired_arg \
    $point_arg \
    $more_args

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
