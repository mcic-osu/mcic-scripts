#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=iqtree
#SBATCH --output=slurm-iqtree-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Construct a phylogenetic tree from a FASTA alignment using IQ-tree"
MODULE=miniconda3
CONDA=/fs/ess/PAS0471/jelmer/conda/iqtree
SCRIPT_VERSION="2023-07-22"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY=iqtree
TOOL_NAME=IQ-tree
TOOL_DOCS=http://www.iqtree.org/doc/
VERSION_COMMAND="$TOOL_BINARY --version"

# Parameter defaults
auto_cores=false                # Don't use IQ-tree's 'AUTO' core mode
ufboot= && boot_arg=            # No bootstrapping by default
model= && model_arg=            # Use IQ-tree's default model (MFP => Pick model)

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
    echo "  - Basic usage:"
    echo "      sbatch $0 -i results/alignment/COI_aligned.fa -p results/iqtree/COI -b 1000"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input FASTA file (alignment) -- should contain multiple, aligned sequences"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --out_prefix    <str>   Output file prefix                          [default: basename of input file]"
    echo "  --model         <str>   Mutation model                              [default: IQ-tree's default = MFP = Pick model]"
    echo "  --ufboot        <int>   Nr of ultrafast bootstraps                  [default: no bootstrapping]"
    echo "  --auto_cores            Use IQ-tree's 'AUTO' core mode              [default: use nr of cores for batch job]"
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
outdir=
out_prefix=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --out_prefix )      shift && out_prefix=$1 ;;
        --model )           shift && model=$1 ;;
        --ufboot )          shift && ufboot=$1 ;;
        --auto_cores )      auto_cores=true ;;
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
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Logging files and dirs
LOG_DIR="$outdir"/logs
VERSION_FILE="$LOG_DIR"/version.txt
CONDA_YML="$LOG_DIR"/conda_env.yml
ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# Define outputs and settings based on script parameters
mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G   # 80% of available memory in GB
[[ "$auto_cores" == true ]] && threads="AUTO"
[[ -n "$ufboot" ]] && boot_arg="--ufboot $ufboot"
[[ -n "$model" ]] && model_arg="-m $model"
[[ -z "$out_prefix" ]] && out_prefix=$(basename "${infile%.*}")

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input FASTA file:                             $infile"
echo "Output dir:                                   $outdir"
echo "Output file prefix:                           $out_prefix"
echo "Use IQ-tree's 'AUTO' core mode:               $auto_cores"
[[ -n "$model" ]] && echo "Model:                                        $model"
[[ -n "$ufboot" ]] && echo "Number of ultrafast bootstraps:               $ufboot"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Run the tool
log_time "Running $TOOL_NAME..."
iqtree \
    -s "$infile" \
    --prefix "$outdir"/"$out_prefix" \
    $model_arg \
    $boot_arg \
    -nt "$threads" -ntmax "$threads" -mem "$mem_gb" \
    -redo \
    $more_args

#? -m MFP  => Model selection with ModelFinder (is the IQ-tree default, too)
#? -redo   => Will overwrite old results
#? --mem 4G  - memory

#TODO - Consider these options:
#-alrt 1000 --wbtl -alninfo

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
