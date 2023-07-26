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
SCRIPT_VERSION="2023-07-25"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=iqtree
TOOL_NAME=IQ-tree
TOOL_DOCS=http://www.iqtree.org/doc
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
conda_path=/fs/ess/PAS0471/jelmer/conda/iqtree
container_path=
container_url=docker://quay.io/biocontainers/iqtree:2.2.2.7--h21ec9f0_2
env=conda                           # 'conda' or 'container'
dl_container=false
container_dir="$HOME/containers"

# Defaults - tool parameters
auto_cores=false                    # Don't use IQ-tree's 'AUTO' core mode
ufboot= && boot_arg=                # No bootstrapping by default
model= && model_arg=                # Use IQ-tree's default model (MFP => Pick model)

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
    echo "  - Basic usage -- include a 100 uf-bootstraps:"
    echo "      sbatch $0 -i results/COI_aligned.fa -o results/iqtree --ufboot 1000"
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
    echo "  --out_prefix        <str>   Output file prefix                          [default: basename of input file]"
    echo "  --model             <str>   Mutation model                              [default: IQ-tree's default = MFP = Pick model]"
    echo "  --ufboot            <int>   Nr of ultrafast bootstraps                  [default: no bootstrapping]"
    echo "  --auto_cores                Use IQ-tree's 'AUTO' core mode              [default: use nr of cores for batch job]"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
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
infile=
outdir=
out_prefix=
opts=
version_only=false
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --out_prefix )      shift && out_prefix=$1 ;;
        --model )           shift && model=$1 ;;
        --ufboot )          shift && ufboot=$1 ;;
        --auto_cores )      auto_cores=true ;;
        --opts )            shift && opts=$1 ;;
        --env )             shift && env=$1 ;;
        --dl_container )    dl_container=true ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 && dl_container=true ;;
        -h | --help )       script_help; exit 0 ;;
        -v )                script_version; exit 0 ;;
        --version )         version_only=true ;;
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
[[ -n "$ufboot" ]] && boot_arg="--ufboot $ufboot"
[[ -n "$model" ]] && model_arg="-m $model"
[[ -z "$out_prefix" ]] && out_prefix=$(basename "${infile%.*}")

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:            $all_opts"
echo "Input FASTA file:                             $infile"
echo "Output dir:                                   $outdir"
echo "Output file prefix:                           $out_prefix"
echo "Use IQ-tree's 'AUTO' core mode:               $auto_cores"
[[ -n "$model" ]] && echo "Model:                                        $model"
[[ -n "$ufboot" ]] && echo "Number of ultrafast bootstraps:               $ufboot"
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
    $model_arg \
    $boot_arg \
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
