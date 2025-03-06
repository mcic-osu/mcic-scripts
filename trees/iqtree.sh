#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=8:00:00
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
SCRIPT_VERSION="2023-10-23"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=iqtree
TOOL_NAME=IQ-tree
TOOL_DOCS=http://www.iqtree.org/doc
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
conda_path=/fs/ess/PAS0471/jelmer/conda/iqtree
container_path=
container_url=docker://quay.io/biocontainers/iqtree:2.2.2.7--h21ec9f0_2
env_type=conda                           # 'conda' or 'container'
dl_container=false
container_dir="$HOME/containers"

# Defaults - tool parameters
auto_cores=false                    # Don't use IQ-tree's 'AUTO' core mode
nboot= && boot_opt=                 # No bootstrapping by default
model= && model_opt=                # Use IQ-tree's default model (MFP => Pick model)
modelset= && modelset_opt=          # Restrict ModelFinder model search to a set of models
fast=false && fast_opt=             # Run IQ-Tree with "--fast" option
root= && root_opt=                  # Pre-specified root/outgroup

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo "                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage -- include 1000 ultrafast-bootstraps:"
    echo "      sbatch $0 -i results/COI_aligned.fa -o results/iqtree --nboot 1000"
    echo "  - Also perform dating on the tree (see http://www.iqtree.org/doc/Dating):"
    echo "      sbatch $0 -i aln.fa -o results/iqtree --opts '--date metadata/dates.tsv --date-ci 100'"
    echo "  - Date a previously inferred tree (see http://www.iqtree.org/doc/Dating):"
    echo "      sbatch $0 -i aln.fa -o results/iqtree --model HKY --opts '--date metadata/dates.tsv --date-ci 100 -te results/iqtree/tree.treefile'"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input FASTA file (alignment) -- should contain multiple, aligned sequences"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --out_prefix        <str>   Output file prefix                      [default: basename of input file]"
    echo "  --root              <str>   Specify outgroup/root sample ID         [default: none]"
    echo "  --fast                      Run IQ-Tree in fast mode                [default: $fast]"    
    echo "  --model             <str>   Mutation model                          [default: IQ-tree's default = MFP => Pick model]"
    echo "  --modelset          <str>   Restrict ModelFinder search to a set of models [default: off]"
    echo "  --nboot             <int>   Nr of ultrafast bootstraps              [default: no bootstrapping]"
    echo "  --auto_cores                Use IQ-tree's 'AUTO' core mode          [default: use nr of cores for batch job]"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env_type               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env_type]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or --dl_container is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: $dl_container]"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
    echo
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
        wget "$FUNCTION_SCRIPT_URL" -O "$function_script"
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
version_only=false
infile=
outdir=
out_prefix=
opts=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --out_prefix )      shift && out_prefix=$1 ;;
        --modelset )        shift && modelset=$1 ;;
        --model )           shift && model=$1 ;;
        --nboot )           shift && nboot=$1 ;;
        --root )            shift && root=$1 ;;
        --auto_cores )      auto_cores=true ;;
        --fast )            fast=true && fast_opt="-fast" ;;
        --opts )            shift && opts=$1 ;;
        --env_type )             shift && env_type=$1 ;;
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
set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G   # 80% of available memory in GB
[[ "$auto_cores" == true ]] && threads="AUTO"
[[ -n "$nboot" ]] && boot_opt="--ufboot $nboot"
[[ -n "$model" ]] && model_opt="-m $model"
[[ -n "$modelset" ]] && modelset_opt="-mset $modelset"
[[ -z "$out_prefix" ]] && out_prefix=$(basename "${infile%.*}")
[[ -n "$root" ]] && root_opt="-o $root"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:            $all_opts"
echo "Input FASTA file:                             $infile"
echo "Output dir:                                   $outdir"
echo "Output file prefix:                           $out_prefix"
echo "Run IQ-Tree in fast mode:                     $fast"
echo "Use IQ-Tree's 'AUTO' core mode:               $auto_cores"
[[ -n "$modelset" ]] && echo "Model set for MFP:                            $modelset"
[[ -n "$model" ]] && echo "Model:                                        $model"
[[ -n "$nboot" ]] && echo "Number of ultrafast bootstraps:               $nboot"
[[ -n "$root" ]] && echo "Root/outgroup:                                $root"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:            $opts"
log_time "Listing the input file(s):"
ls -lh "$infile"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    -s "$infile" \
    --prefix "$outdir"/"$out_prefix" \
    $model_opt \
    $modelset_opt \
    $fast_opt \
    $boot_opt \
    $root_opt \
    -nt "$threads" -ntmax "$threads" -mem "$mem_gb" \
    -redo \
    $opts

#? Options used:
#> -m MFP       => Model selection with ModelFinder (is the IQ-tree default, too)
#> -redo        => Will overwrite old results

#? Consider these options, too:
#> -alrt 1000 --wbtl -alninfo

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
