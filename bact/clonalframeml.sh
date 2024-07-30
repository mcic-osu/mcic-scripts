#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=clonalframeml
#SBATCH --output=slurm-clonalframeml-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run ClonalFrameML to infer recombination (HGT) in bacterial genomes,
and run maskrc-svg to mask recombinant regions"
MODULE=miniconda3
CONDA=/fs/ess/PAS0471/jelmer/conda/clonalframeml
CONTAINER=/fs/ess/PAS0471/containers/depot.galaxyproject.org-singularity-mulled-v2-f5c68f1508671d5744655da9b0e8b609098f4138-7e089189af7822a6a18245830639dbfe11a4c277-0.img
SCRIPT_VERSION="2023-07-22"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY=ClonalFrameML
TOOL_NAME=ClonalFrameML
TOOL_DOCS=https://github.com/xavierdidelot/ClonalFrameML
VERSION_COMMAND="$TOOL_BINARY | head -n 1"

# Defaults - parameters
tree_tool=iqtree            #! NOTE: ClonalFrameML may not accept FastTree trees?! Have had problems with this at least
force_tree=false            # If tree is found, don't rerun tree-building

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
    echo "      sbatch $0 -i results/roary/core_gene_alignment.aln -o results/clonalframeml"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile     <file>  Input alignment file from Roary/Panaroo or similar"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --out_prefix    <str>   Output file prefix                          [default: basename of input file]"
    echo "  --tree_tool     <str>   Tool to build initial tree, 'iqtree' or 'fasttree' [default: 'iqtree']"
    echo "  --force_tree            Even if the initial tree file is found, rebuild it [default: use existing tree]"
    echo "  --more_args     <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - ClonalFrameML: $TOOL_DOCS"
    echo "  - maskrc-svg: https://github.com/kwongj/maskrc-svg"
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
        --tree_tool )       shift && tree_tool=$1 ;;
        --force_tree )      force_tree=true ;;
        --out_prefix )      shift && out_prefix=$1 ;;
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

# Define settings and outputs based on script parameters
mem_gb=4G && [[ "$IS_SLURM" == true ]] && mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G
[[ -z "$out_prefix" ]] && out_prefix=$(basename "${infile%.*}")

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options/arguments passed to this script:  $all_args"
echo "Input alignment file:                         $infile"
echo "Output dir:                                   $outdir"
echo "Output file prefix:                           $out_prefix"
echo "Tool to build initial tree:                   $tree_tool"
echo "Force tree-building even if file exists:      $force_tree"
[[ -n $more_args ]] && echo "Additional arguments for $TOOL_NAME:          $more_args"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Get starting tree
if [[ ! -f "$outdir"/"$out_prefix".treefile || "$force_tree" == true ]]; then
    log_time "Running $tree_tool to get a starting tree for ClonalFrameML..."
    
    if [[ "$tree_tool" == "iqtree" ]]; then
        runstats iqtree \
            -s "$infile" --prefix "$outdir"/"$out_prefix" \
            -m MFP -fast -redo -nt "$threads" -ntmax "$threads" -mem "$mem_gb"
    elif [[ "$tree_tool" == "fasttree" ]]; then
        runstats fasttree -nt -gtr < "$infile" > "$outdir"/"$out_prefix".treefile
    fi
else
    log_time "Tree file $outdir/$out_prefix.treefile exists, skipping tree-building..."
fi
log_time "Tree file:"
ls -lh "$outdir"/"$out_prefix".treefile

# Run ClonalFrameML
if [[ ! -s "$outdir"/"$out_prefix".labelled_tree.newick ]]; then
    log_time "Running ClonalFrameML..."
    runstats $TOOL_BINARY \
        "$outdir"/"$out_prefix".treefile \
        "$infile" \
        "$outdir"/"$out_prefix" \
        $more_args
else
    log_time "ClonalFrameML output file $outdir/$out_prefix.labelled_tree.newick exists, skipping ClonalFrameML..."
fi

# Run maskrc-svg.py to HGT-mask the alignment
# Have to run maskrc-svg.py using a container because the Conda env doesn't work
log_time "Running maskrc-svg.py..."
runstats singularity exec "$CONTAINER" maskrc-svg.py \
    --aln "$infile" \
    --symbol '-' \
    --out "$outdir"/"$out_prefix".masked.aln \
    "$outdir"/"$out_prefix"

#? --symbol CHAR        symbol to use for masking (default="?")
#? The positional argument is the prefix used for ClonalFrameML (or Gubbins) output files

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$VERSION_COMMAND" "$VERSION_FILE" "$ENV_FILE" "$IS_SLURM" \
    "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL"
