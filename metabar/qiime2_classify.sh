#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=25
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=qiime-classify
#SBATCH --output=slurm-qiime-classify-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Classify ASV sequences with Qiime2's feature-classifier classify-sklearn"
SCRIPT_VERSION="2025-02-26"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=qiime
TOOL_NAME=Qiime
TOOL_DOCS=https://docs.qiime2.org/
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/qiime2-amplicon-2024.10
container_url=
container_dir="$HOME/containers"
version_only=false                 # When true, just print tool & script version info and exit

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
    echo "  - Basic usage example:"
    echo "      sbatch $0 --query_fasta ASVf.fa --ref_fasta REF.fa --ref_tax TAX.tsv -o results/qiime --prefix eukaryome"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --query_fasta       <file>  Input query (ASV) FASTA file"
    echo "  --ref_fasta         <file>  Input reference FASTA file"
    echo "  --ref_tax           <file>  Input reference TSV taxonomy file"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  --prefix            <file>  Output file prefix"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more_opts         <str>   Quoted string with additional options for"
    echo "                              the feature-classifier classify-sklearn command"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --container_path    <file>  Pre-existing Singularity container image file (.sif) to use"
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
    function_script_name="$(basename "$FUNCTION_SCRIPT_URL")"
    function_script="$script_dir"/../dev/"$function_script_name"

    # Download the function script if needed, then source it
    if [[ -f "$function_script" ]]; then
        source "$function_script"
    elif [[ ! -f "$function_script_name" ]]; then
        echo "Can't find script with Bash functions ($function_script_name), downloading from GitHub..."
        wget -q "$FUNCTION_SCRIPT_URL" -O "$function_script_name"
        source "$function_script_name"
    else
        source "$function_script_name"
    fi
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
query_fasta=
ref_fasta=
ref_tax=
prefix=
outdir=
more_opts=
threads=
container_path=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --query_fasta )     shift && query_fasta=$1 ;;
        --ref_fasta )       shift && ref_fasta=$1 ;;
        --ref_tax )         shift && ref_tax=$1 ;;
        --prefix )          shift && prefix=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env )             shift && env=$1 ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 ;;
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
load_env "$conda_path" "$container_path"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$query_fasta" ]] && die "No input query FASTA file specified, do so with --query_fasta" "$all_opts"
[[ -z "$ref_fasta" ]] && die "No input reference FASTA file specified, do so with --ref_fasta" "$all_opts"
[[ -z "$ref_tax" ]] && die "No input reference taxonomy file specified, do so with --ref_tax" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ -z "$prefix" ]] && die "No output file prefix specified, do so with --prefix" "$all_opts"
[[ ! -f "$query_fasta" ]] && die "Input query FASTA file $query_fasta does not exist"
[[ ! -f "$ref_fasta" ]] && die "Input reference FASTA file $ref_fasta does not exist"
[[ ! -f "$ref_tax" ]] && die "Input reference taxonomy file $ref_tax does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Working directory:                        $PWD"
echo
echo "Input query FASTA file:                   $query_fasta"
echo "Input reference FASTA file:               $ref_fasta"
echo "Input reference taxonomy file:            $ref_tax"
echo "Output dir:                               $outdir"
echo "Output prefix:                            $prefix"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$query_fasta" "$ref_fasta" "$ref_tax"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Importing the ASV sequences
log_time "Importing the ASV sequences..."
runstats $CONTAINER_PREFIX $TOOL_BINARY tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$query_fasta" \
    --output-path "$outdir"/"$prefix"_query-seq.qza

# Importing the reference taxonomy
log_time "Importing the reference taxonomy..."
runstats $CONTAINER_PREFIX $TOOL_BINARY tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-path "$ref_tax" \
    --output-path "$outdir"/"$prefix"_ref-tax.qza

# Importing the reference sequence
log_time "Importing the reference sequence..."
runstats $CONTAINER_PREFIX $TOOL_BINARY tools import \
    --type 'FeatureData[Sequence]' \
    --input-path "$ref_fasta" \
    --output-path "$outdir"/"$prefix"_ref-seq.qza

# Train the classifier
log_time "Training the classifier..."
runstats $CONTAINER_PREFIX $TOOL_BINARY feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$outdir"/"$prefix"_ref-seq.qza \
    --i-reference-taxonomy "$outdir"/"$prefix"_ref-tax.qza \
    --o-classifier "$outdir"/"$prefix"_classifier.qza

# Run the classifier
log_time "Running the classifier..."
runstats $CONTAINER_PREFIX $TOOL_BINARY feature-classifier classify-sklearn \
    --i-classifier "$outdir"/"$prefix"_classifier.qza \
    --i-reads "$outdir"/"$prefix"_query-seq.qza \
    --o-classification "$outdir"/"$prefix"_taxonomy.qza \
    --p-n-jobs "$threads" \
    $more_opts

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
