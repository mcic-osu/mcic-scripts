#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=diamond
#SBATCH --output=slurm-diamond-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run DIAMOND to perform fast BLAST-like alignment of proteins"
SCRIPT_VERSION="2023-12-09"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=diamond
TOOL_NAME=Diamond
TOOL_DOCS=https://github.com/bbuchfink/diamond/wiki
VERSION_COMMAND="$TOOL_BINARY version"

# Defaults - generics
env_type=container                       # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/diamond
container_path=
container_url=docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0
dl_container=false
container_dir="$HOME/containers"
strict_bash=true
version_only=false                 # When true, just print tool & script version info and exit

# Defaults - tool parameters
blast_type=blastp                   # Or blastx
out_format="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp slen stitle"
max_target_seqs=25                  # Same as DIAMOND default
evalue="0.001"                      # E-value threshold
pct_id=0                            # % identity threshold (empty => no threshold)
pct_qcov=0                          # Threshold for % of query covered by the alignment length
pct_scov=0                          # Threshold for % of subject covered by the alignment length
add_header=true                     # Add column header to final BLAST output file

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
    echo "      sbatch $0 -i query.fa -o results/diamond --db data/diamond_db.dmnd"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input FASTA file (can contain one or more sequences)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  --db                <str>   Diamond DB '.dmnd' file (can create this with diamond_db.sh)"
    echo
    echo "GENERAL OPTIONS (OPTIONAL):"
    echo "  --blast_type        <str>   BLAST type: 'blastp' or 'blastx'        [default: $blast_type]"
    echo "  --sens              <str>   Sensitivity: one of 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'"
    echo "  --out_format        <str>   Output format string. NOTE: changing this may mess up output filtering steps, which rely on the default format"
    echo "                                [default: $out_format]"
    echo "  --no_header                 Don't add column headers to final output TSV file [default: add]"
    echo "                                The header won't be added to the raw output file, which can be used for filtering"
    echo "  --more_opts         <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "THRESHOLD AND FILTERING OPTIONS (OPTIONAL):"
    echo "  --max_target_seqs   <int>   Max. nr of target sequences to keep                     [default: DIAMOND default (=25)]"
    echo "  --evalue            <num>   E-value threshold in scientific notation                [default: $evalue]"
    echo "  --pct_id            <int>   Percentage identity threshold                           [default: $pct_id]"
    echo "  --pct_qcov          <int>   Threshold for % of query covered by the alignment       [default: $pct_qcov]"
    echo "  --pct_qcov          <int>   Threshold for % of query covered by the alignment       [default: $pct_scov]"
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
db=
sensitivity=
header_opt=
more_opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --sens )            shift && sensitivity=$1 ;;
        --out_format )      shift && out_format=$1 ;;
        --no_header )       add_header=false ;;
        --max_target_seqs ) shift && max_target_seqs=$1 ;;
        --db )              shift && db=$1 ;;
        --blast_type )      shift && blast_type=$1 ;;
        --evalue )          shift && evalue=$1 ;;
        --pct_id )          shift && pct_id=$1 ;;
        --pct_qcov )         shift && pct_qcov=$1 ;;
        --pct_scov )         shift && pct_scov=$1 ;;
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
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$db" ]] && die "No database file specified, do so with --db" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -f "$db" ]] && die "Database file $db does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
outfile="$outdir"/diamond_out.tsv
[[ "$add_header" == true ]] && header_opt="--header"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input file:                               $infile"
echo "Output dir:                               $outdir"
echo "DIAMOND db:                               $db"
echo
echo "BLAST type:                               $blast_type"
echo "Add column header to output?              $add_header"
[[ -n "$sensitivity" ]] && echo "Sensitivity:                              $sensitivity"
echo "Evalue threshold:                         $evalue"
[[ -n "$pct_id" ]] && echo "Percent identity threshold:               $pct_id"
[[ -n "$pct_qcov" ]] && echo "Percent query coverage threshold:         $pct_qcov"
[[ -n "$pct_scov" ]] && echo "Percent subject coverage threshold:       $pct_scov"
[[ -n "$max_target_seqs" ]] && echo "Max. nr. of target sequences:             $max_target_seqs"
echo "Number of queries in the input file:      $(grep -c "^>" "$infile")"
log_time "Listing the input file(s):"
ls -lh "$infile"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY $blast_type \
    --db "$db" \
    --query "$infile" \
    --out "$outfile" \
    --outfmt $out_format \
    --max-target-seqs "$max_target_seqs" \
    --evalue "$evalue" \
    --id "$pct_id" \
    --query-cover "$pct_qcov" \
    --subject-cover "$pct_scov" \
    --"$sensitivity" \
    --threads "$threads" \
    $header_opt \
    $more_opts

# ==============================================================================
#                           REPORT & WRAP UP
# ==============================================================================
# Report some basic stats on the output
if [[ "$add_header" == false ]]; then
    n_hits=$(wc -l < "$outfile")
    n_queries=$(cut -f 1 "$outfile" | sort -u | wc -l)
    n_subjects=$(cut -f 2 "$outfile" | sort -u | wc -l)
else
    n_hits=$(tail -n+4 "$outfile" | wc -l)
    n_queries=$(tail -n+4 "$outfile" | sort -u | wc -l)
    n_subjects=$(tail -n+4 "$outfile" | cut -f 2 | sort -u | wc -l)
fi
log_time "Number of hits in the final output file: $n_hits"
log_time "Number of distinct queries in the final output file: $n_queries"
log_time "Number of distinct subjects in the final output file: $n_subjects"

# Final logging
log_time "Listing the output file:"
ls -lh "$outfile"
final_reporting "$LOG_DIR"
