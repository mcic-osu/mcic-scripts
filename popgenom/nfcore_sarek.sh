#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=nfc_sarek
#SBATCH --output=slurm-nfc_sarek-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run the Nextflow/nf-core Sarek pipeline for non-model organism genomic variant callling"
SCRIPT_VERSION="2025-01-18"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY="nextflow run"
export TOOL_NAME="nextflow"
VERSION_COMMAND="nextflow -v"

# Constants - specific
WORKFLOW_NAME=nf-core/sarek                            # The name of the nf-core workflow
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config # Nextflow <=> OSC config file
NONMODEL_OPTS="--igenomes_ignore --skip_tools baserecalibrator"

# Defaults - pipeline parameters
workflow_version=3.5.0                                 # The version of the nf-core workflow
resume=true && resume_arg="-resume"                    # Resume the workflow from wherever it left off

# Defaults - infrastructure
conda_path=/fs/project/PAS0471/jelmer/conda/nextflow   # Conda environment with Nextflow & nf-core tools
osc_account=PAS0471                                    # If the script is submitted with another project, this will be updated (line below)
[[ -n $SLURM_JOB_ACCOUNT ]] && osc_account=$(echo "$SLURM_JOB_ACCOUNT" | tr "[:lower:]" "[:upper:]")
container_dir=/fs/scratch/"$osc_account"/containers    # The workflow will download containers to this dir
work_dir=/fs/scratch/"$osc_account"/$USER/nfc-sarek    # 'work dir' for initial outputs (selected, final outputs go to the outdir)
profile="singularity"                                  # 'singularity' to have the workflow use containers (alternatively, 'conda')
version_only=false                                     # When true, just print tool & script version info and exit 

# ==============================================================================
#                                FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLES:"
    echo "  - Basic usage (always submit your scripts to SLURM with 'sbatch'):"
    echo "      sbatch $0 --samplesheet metadata/samplesheet.csv --genome data/assembly.fna -o results/sarek --tools freebayes"
    echo "  - Use the Freebayes variant caller in addition to mpileup:"
    echo "      sbatch $0 --samplesheet metadata/samplesheet.csv --genome data/assembly.fna -o results/sarek --tools freebayes,mpileup"
    echo "  - Example of using '--more_args': save the BAM files"
    echo "      sbatch $0 --samplesheet metadata/samplesheet.csv --genome data/assembly.fna -o results/sarek --more_args '--save_mapped --save_output_as_bam'"
    echo "  - To just print the help message for this script:"
    echo "      bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--samplesheet    <file>  Input samplesheet CSV file (see https://nf-co.re/sarek/usage)"
    echo "  --genome            <file>  Reference genome FASTA file"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --tools             <str>   Comma-separated list of one or more tools for variant calling and annotation"
    echo "                              See https://nf-co.re/sarek/parameters#tools for options   [default: 'mpileup']"
    echo "  --workflow_version  <str>   Nf-core rnaseq workflow version         [default: $workflow_version]"
    echo "  --more_opts         <str>   Additional workflow parameters, for usage see example above"
    echo "                              Available workflow parameters are listed at https://nf-co.re/sarek/parameters"
    echo
    echo "NEXTFLOW-RELATED OPTIONS:"
    echo "  --restart                   Restart workflow from the beginning     [default: resume workflow if possible]"
    echo "  --container_dir     <dir>   Directory with container images         [default: $container_dir]"
    echo "                                - Required container images will be downloaded here when not already present" 
    echo "  --config            <file>  Additional Nextflow config file         [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --profile            <str>  Name of a 'profile' from a config file that should be used [default: $profile]"
    echo "  --work_dir           <dir>  Scratch (work) dir for the workflow     [default: $work_dir]"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
    echo "  --version                   Print the version of Nextflow and exit"
    echo
    echo "UTILITY OPTIONS"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo
    echo "HARDCODED WORKFLOW OPTIONS & DEFAULTS THAT DIFFER FROM THE WORKFLOW DEFAULTS:"
    echo "  - The option '--aligner star_salmon' is hardcoded"
    echo
    echo "NF-CORE SAREK WORKFLOW DOCUMENTATION:"
    echo "   - https://nf-co.re/sarek"
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
    function_script=$(realpath "$script_dir"/../dev/bash_functions2.sh)
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions2.sh
    fi
    source "$function_script"
}

nextflow_setup() {
    # Singularity container dir - any downloaded containers will be stored here;
    # if the required container is already there, it won't be re-downloaded
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
    export NXF_SINGULARITY_LIBRARYDIR="$container_dir"
    mkdir -p "$NXF_SINGULARITY_CACHEDIR"

    # Limit memory for Nextflow main process - see https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html
    export NXF_OPTS='-Xms1g -Xmx4g'
}

# ==============================================================================
#                          INFRASTRUCTURE SETUP I
# ==============================================================================
# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script $IS_SLURM

# ==============================================================================
#                     PARSE COMMAND-LINE OPTIONS
# ==============================================================================
samplesheet=
ref_fasta=
outdir=
config_file=
more_opts=
threads=

# Parse command-line options
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )             shift && outdir=$1 ;;
        --samplesheet )             shift && samplesheet=$1 ;;
        --ref_fasta )               shift && ref_fasta=$1 ;;
        --more_opts )               shift && more_opts=$1 ;;
        --workflow_version )        shift && workflow_version=$1 ;;
        --container_dir )           shift && container_dir=$1 ;;
        --config | -config )        shift && config_file=$1 ;;
        --profile | -profile )      shift && profile=$1 ;;
        --work_dir | -work-dir )    shift && work_dir=$1 ;;
        --restart | -restart )      resume=false && resume_arg= ;;
        -h | --help )               script_help; exit 0 ;;
        -v )                        script_version; exit 0 ;;
        --version )                 version_only=true ;;
        * )                         script_help; die "Invalid option $1";;
    esac
    shift
done

# Check input
[[ -z "$samplesheet" ]] && die "Please specify a samplesheet with --samplesheet" "$all_opts"
[[ -z "$ref_fasta" ]] && die "Please specify a reference genome FASTA file with --ref_fasta" "$all_opts"
[[ -z "$outdir" ]] && die "Please specify an output dir with -o/--outdir" "$all_opts"
[[ ! -f "$samplesheet" ]] && die "Samplesheet $samplesheet does not exist"
[[ ! -f "$ref_fasta" ]] && die "Reference FASTA file $ref_fasta does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Bash strict settings
set -euo pipefail

# Load software
load_env "$conda_path"
nextflow_setup
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Build the config argument
OSC_CONFIG="$outdir"/$(basename "$OSC_CONFIG_URL")
config_arg="-c $OSC_CONFIG"
[[ -n "$config_file" ]] && config_arg="$config_arg -c ${config_file/,/ -c }"

# Other output dirs
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
trace_dir="$outdir"/pipeline_info

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script: $all_opts"
echo
echo "INPUT AND OUTPUT:"
echo "Sample sheet:                     $samplesheet"
echo "Reference genome FASTA file:      $ref_fasta"
echo "Output dir:                       $outdir"
[[ -n "$more_opts" ]] && echo "Additional options:               $more_opts"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run?              $resume"
echo "Container dir:                    $container_dir"
echo "Scratch (work) dir:               $work_dir"
echo "Config 'profile':                 $profile"
[[ -n "$config_file" ]] && echo "Additional config file:           $config_file"
echo "Config file argument:             $config_arg"
echo "=========================================================================="
echo "Listing the input files and showing the first lines of the samplesheet:"
ls -lh "$ref_fasta"
echo
head "$samplesheet"
echo "=========================================================================="
set_threads "$IS_SLURM"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                              RUN
# ==============================================================================
# Make necessary dirs
log_time "Creating the output directories..."
mkdir -pv "$work_dir" "$container_dir" "$outdir"/logs "$trace_dir"

# Download the OSC config file
if [[ ! -f "$OSC_CONFIG" ]]; then
    log_time "Downloading the mcic-scripts Nextflow OSC config file to $OSC_CONFIG..."
    wget -q -O "$OSC_CONFIG" "$OSC_CONFIG_URL"
fi

# Modify the config file so it has the correct OSC project/account
if [[ "$osc_account" != "PAS0471" ]]; then
    sed -i "s/--account=PAS0471/--account=$osc_account/" "$OSC_CONFIG"
fi

# Run the workflow
log_time "Starting the workflow.."
runstats $TOOL_BINARY $WORKFLOW_NAME \
    -r $workflow_version \
    --input "$samplesheet" \
    --outdir "$outdir" \
    --fasta "$ref_fasta" \
    --tools "$tools" \
    $NONMODEL_OPTS \
    -work-dir "$work_dir" \
    -with-report "$trace_dir"/report.html \
    -with-trace "$trace_dir"/trace.txt \
    -with-timeline "$trace_dir"/timeline.html \
    -with-dag "$trace_dir"/dag.png \
    -ansi-log false \
    -profile "$profile" \
    $config_arg \
    $resume_arg \
    $more_opts

# Report
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
