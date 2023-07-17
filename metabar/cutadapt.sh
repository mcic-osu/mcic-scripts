#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=cutadapt
#SBATCH --output=slurm-cutadapt-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Run Cutadapt to remove metabarcoding primers for a single pair of FASTQ files
The script will compute and use the reverse complements of all primers as well."
readonly MODULE=miniconda3
readonly CONDA=/fs/ess/PAS0471/jelmer/conda/cutadapt
readonly SCRIPT_VERSION="2023-07-17"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY=cutadapt
readonly TOOL_NAME=CutAdapt
readonly TOOL_DOCS=https://cutadapt.readthedocs.io/en/stable
readonly VERSION_COMMAND="$TOOL_BINARY --version"

# Parameter defaults
discard_untrimmed=true

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage (always submit your scripts to SLURM with 'sbatch'):"
    echo "      sbatch $0 -i data/sample1_R1.fastq.gz -o results/cutadapt -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT"
    echo "  - Using a primer file instead:"
    echo "      sbatch $0 -i data/sample1_R1.fastq.gz -o results/cutadapt --primer_file metadata/primers.txt"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1             <file>  Input R1 FASTQ file (name of R2 will be inferred)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "There are two ways of specifying primers: (1) with '-f' and '-r' or (2) with '--primer_file' (use the latter if you have multiple primer pairs):"
    echo "  -f/--primer_f       <str>   Forward primer sequence (use in combination with -r)"
    echo "  -r/--primer_r       <str>   Reverse primer sequence (use in combination with -f)"
    echo "  --primer_file       <file>  File with primer sequences, one pair per line separated by a space (*alternative* to using -f and -r)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --keep_untrimmed            Don't discard untrimmed sequences (i.e. those with no primers) [default: discard]"
    echo "  --more_args         <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "HARDCODED OPTIONS:"
    echo "  - The CutAdapt option '--pair-filter=any' is always used."
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
R1=
outdir=
discard_arg=
primer_f=
primer_r=
primer_file=
primer_arg=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && readonly R1=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        -f | --primer_f )   shift && primer_f=$1 ;;
        -r | --primer_r )   shift && primer_r=$1 ;;
        --primer_file )     shift && primer_file=$1 ;;
        --keep_untrimmed )  discard_untrimmed=false ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -v )                script_version; exit 0 ;;
        -h | --help )       script_help; exit 0 ;;
        --version )         load_env "$MODULE" "$CONDA"
                            tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$R1" ]] && die "No input file specified, do so with -i/--R1" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$R1" ]] && die "Input file $R1 does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Logging files and dirs
readonly LOG_DIR="$outdir"/logs
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# Define outputs based on script parameters
[[ "$discard_untrimmed" == true ]] && discard_arg="--discard-untrimmed"

## Determine input dir and R2 file
indir=$(dirname "$R1")
R1_base=$(basename "$R1")
file_ext=$(echo "$R1_base" | sed -E 's/.*(.fastq|.fq|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2="$indir"/${R1_base/$R1_suffix/$R2_suffix}
R1_basename=$(basename "$R1")
R2_basename=$(basename "$R2")
[[ ! -f "$R2" ]] && die "Input FASTQ file $R2 not found"
[[ "$indir" == "$outdir" ]] && die "Input dir should not be the same as output dir ($indir)"
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input R1 FASTQ file:                          $R1"
echo "Input R2 FASTQ file:                          $R2"
echo "Output dir:                                   $outdir"
echo "Discard untrimmed (-d):                       $discard_untrimmed"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$R1" "$R2"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               DEFINE PRIMERS
# ==============================================================================
if [[ -z "$primer_file" ]]; then
    log_time "Using the forward and reverse primers provided as arguments to the script..."
    [[ -z "$primer_f" ]] && die "No forward primer (-f) provided"
    [[ -z "$primer_r" ]] && die "No reverse primer (-r) provided"

    primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
    primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

    primer_arg="-a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"

    echo "Forward primer (-f):              $primer_f"
    echo "Reverse primer (-r):              $primer_r"
    echo "Forward primer - rev. comp.:      $primer_f_rc"
    echo "Reverse primer - rev. comp.:      $primer_r_rc"

else
    log_time "Using primer file $primer_file to read primers..."
    [[ ! -f "$primer_file" ]] && die "Primer file $primer_file not found"

    while read -r primer_f primer_r; do

        primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
        primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

        echo "Forward primer (-f):          $primer_f"
        echo "Reverse primer (-r):          $primer_r"
        echo "Forward primer - rev. comp.:  $primer_f_rc"
        echo "Reverse primer - rev. comp.:  $primer_r_rc"

        primer_arg="$primer_arg -a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"
        primer_arg=$(echo "$primer_arg" | sed -E 's/^ +//') # Remove leading whitespace

    done <"$primer_file"
fi
echo "Primer argument:                  $primer_arg"

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
        $primer_arg \
        --output "$outdir"/"$R1_basename" \
        --paired-output "$outdir"/"$R2_basename" \
        --pair-filter=any \
        $discard_arg \
        --cores "$threads" \
        $more_args \
        "$R1" "$R2"

#? --pair-filter=any: Remove pair if one read is filtered (=Default)

# List the output, report version, etc
log_time "Listing the output files:"
ls -lhd "$(realpath "$outdir")"/"$sample_id"*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
