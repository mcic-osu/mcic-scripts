#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=fastp
#SBATCH --output=slurm-fastp-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Run fastp to preprocess FASTQ files"
readonly MODULE=miniconda3
readonly CONDA=/fs/ess/PAS0471/jelmer/conda/fastp
readonly SCRIPT_VERSION="2023-07-15"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY=fastp
readonly TOOL_NAME=fastp
readonly TOOL_DOCS=https://github.com/OpenGene/fastp
readonly TOOL_PAPER=https://academic.oup.com/bioinformatics/article/34/17/i884/5093234
readonly VERSION_COMMAND="$TOOL_BINARY --version"

# Parameter defaults
save_unpaired=false
single_end=false

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
    echo "      sbatch $0 -i data/A_R1.fastq.gz -o results/fastp"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  Input (R1/forward) FASTQ file (name of R2 file will be inferred)"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo "                            Output FASTQ files will have the same names as the input files"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --single_end            Sequences are single-end: don't look for R2 file"
    echo "  --save_unpaired         Save all unpaired (orphaned due to QC) sequences in a single output file"
    echo "                            This FASTQ file's name will end in '_unpaired', e.g. 'sampleA_unpaired.fastq.gz"
    echo "  --more_args     <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Online docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
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
R1_in=
R2_in=
outdir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && readonly R1_in=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        --single_end )      readonly single_end=true ;;
        --save_unpaired )   readonly save_unpaired=true ;;
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
[[ -z "$R1_in" ]] && die "No input file specified, do so with -i/--R1" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$R1_in" ]] && die "Input file $R1_in does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Logging files and dirs
readonly LOG_DIR="$outdir"/logs
readonly REPORT_DIR="$outdir"/reports
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR" "$REPORT_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# Define outputs based on script parameters
file_ext=$(basename "$R1_in" | sed -E 's/.*(.fasta|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1_in" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
sample_id=$(basename "$R1_in" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")
R1_out="$outdir"/"$sample_id""$R1_suffix""$file_ext"

if [[ "$single_end" == false ]]; then
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    [[ ! -f "$R2_in" ]] && die "Input file $R2_in does not exist"
    R2_out="$outdir"/"$sample_id""$R2_suffix""$file_ext"
    R2_arg="--in2 $R2_in --out2 $R2_out"
fi

if [[ "$save_unpaired" == true ]]; then
    [[ "$single_end" == true ]] && die "Can't save unpaired output when input is single-end"
    unpaired_out="$outdir"/"$sample_id"_unpaired"$file_ext"
    unpaired_arg="--unpaired1 $unpaired_out --unpaired2 $unpaired_out"
fi

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input (R1) file:                              $R1_in"
[[ "$single_end" == false ]] && echo "Input R2 file:                                $R2_in"
echo "Output (R1) file:                             $R1_out"
[[ "$single_end" == false ]] && echo "Output R2 file:                               $R2_out"
[[ "$save_unpaired" == true ]] && echo "Output unpaired file:                         $unpaired_out"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$R1_in"
[[ "$single_end" == false ]] && ls -lh "$R2_in"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --in1 "$R1_in" \
    --out1 "$R1_out" \
    $R2_arg \
    $unpaired_arg \
    --json "$REPORT_DIR"/"$sample_id".json \
    --html "$REPORT_DIR"/"$sample_id".html \
    --thread "$threads" \
    $more_args

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"

#? --detect_adapter_for_pe          by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.
