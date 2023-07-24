#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --mem=172G
#SBATCH --cpus-per-task=40
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=pilon
#SBATCH --output=slurm-pilon-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Pilon to polish a genome assembly with Illumina reads"
MODULE=miniconda3
CONDA=/fs/ess/PAS0471/jelmer/conda/pilon-1.24
SCRIPT_VERSION="2023-07-23"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
PILON_JAR=/fs/ess/PAS0471/jelmer/conda/pilon-1.24/share/pilon-1.24-0/pilon.jar
TOOL_NAME=Pilon
TOOL_DOCS=https://github.com/broadinstitute/pilon/wiki
VERSION_COMMAND=

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "          $0 (v. $SCRIPT_VERSION)"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage:"
    echo "      sbatch $0 -i results/assembly.fna --bam_dir results/bwa -o results/pilon"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--assembly   <file>  Genome assembly FASTA file"
    echo "  --bam_dir       <dir>   Directory with BAM files of Illumina reads mapped to the assembly"
    echo "  --outfile       <file>  Output assembly file (dir will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --fix           <str>   What to fix: 'snps'/'indels'/'gaps'/'local'/'all'/'bases'   [default: Pilon default => 'all']"
    echo "                          See the Pilon documentation for details"
    echo "  --more_args     <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
    echo
}

# Function to source the script with Bash functions
source_function_script() {
    local is_slurm=$1

    # Determine the location of this script, and based on that, the function script
    if [[ "$is_slurm" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script=$(realpath "$script_dir"/../dev/bash_functions.sh)
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions.sh
    fi
    source "$function_script"
}

# ==============================================================================
#                          INFRASTRUCTURE SETUP I
# ==============================================================================
# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script $IS_SLURM

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
infile=
outfile=
bam_arg=
fix= && fix_arg=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outfile )    shift && outfile=$1 ;;
        --bam_dir )         shift && bam_dir=$1 ;;
        --fix )             shift && fix=$1 && fix_arg="--fix $fix" ;;
        --more_args )       shift && more_args=$1 ;;
        -v )                script_version; exit 0 ;;
        -h | --help )       script_help; exit 0 ;;
        --version )         load_env "$MODULE" "$CONDA"
                            tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_args"
[[ -z "$outfile" ]] && die "No output file specified, do so with -o/--outfile" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Determine the output dir
outdir=$(dirname "$outfile")

# Logging files and dirs
LOG_DIR="$outdir"/logs
VERSION_FILE="$LOG_DIR"/version.txt
CONDA_YML="$LOG_DIR"/conda_env.yml
ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# Other settings based on script options
for bam in "$bam_dir"/*bam; do bam_arg="$bam_arg --frags $bam"; done
mem=$(( SLURM_MEM_PER_NODE / 1000 ))G # Get amount of RAM in right format
file_ext="${outfile##*.}"
out_prefix=$(basename "$outfile" ."$file_ext")

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input assembly FASTA:                         $infile"
echo "Input BAM dir:                                $bam_dir"
echo "Output assembly FASTA:                        $outfile"
[[ -n "$fix" ]] && echo "What to fix (--fix):                          $fix"
echo "Memory for Java:                              $mem"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input assembly:"
ls -lh "$infile"
log_time "Listing the input BAM files:"
ls -lh "$bam_dir"/*bam
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
runstats java -Xmx"$mem" -jar "$PILON_JAR" \
    --genome "$infile" \
    $bam_arg \
    --outdir "$outdir" \
    --output "$out_prefix" \
    $fix_arg \
    $more_args

#? No --threads arg: running with v 1.24, got: "--threads argument no longer supported; ignoring!"

# Rename the output file if needed
if [[ "$outfile" != "$outdir"/"$out_prefix".fasta ]]; then
    log_time "Renaming the output file..."
    mv -v "$outdir"/"$out_prefix".fasta "$outfile"
fi

# List the output, report version, etc
log_time "Listing the output file:"
ls -lh "$outfile"
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
