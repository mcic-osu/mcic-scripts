#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=sarek
#SBATCH --output=slurm-sarek-%j.out

#TODO - No way of providing annotation via GFF/GTF?

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Run the Nextflow/nf-core Sarek pipeline for non-model organism genomic variant callling"
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA=/fs/project/PAS0471/jelmer/conda/nextflow
readonly SCRIPT_VERSION="1.1"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY="nextflow run"
readonly TOOL_NAME="nf-core Sarek"
readonly TOOL_DOCS=https://nf-co.re/sarek
readonly TOOL_PAPER=http://dx.doi.org/10.12688/f1000research.16665.2

# Constants - specific
#? By default, the workflow will download human 'igenome' data - skip this
#? Then also skip the base score recalibration in the absence of reference data
readonly NONMODEL_OPTS="--igenomes_ignore --skip_tools baserecalibrator"

# Option defaults - Nextflow
workflow_version="3.2.0"  # Sarek version #TODO Update to latest version once container download issue is resolved
workflow_dir=workflows/nfcore-sarek
container_dir=/fs/project/PAS0471/containers
work_dir=/fs/scratch/PAS0471/$USER/nfcore-sarek
profile="singularity"
resume=true && resume_arg="-resume"
download_wf=false
variant_tools=mpileup

# URL to OSC Nextflow config file
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
osc_config=mcic-scripts/nextflow/osc.config  # Will be downloaded if not present here

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
    echo "      sbatch $0 --samplesheet metadata/samplesheet.csv --genome data/assembly.fna -o results/sarek --tools freebayes"
    echo "  - Use the Freebayes variant caller in addition to mpileup:"
    echo "      sbatch $0 --samplesheet metadata/samplesheet.csv --genome data/assembly.fna -o results/sarek --tools freebayes,mpileup"
    echo "  - Example of using '--more_args': save the BAM files"
    echo "      sbatch $0 --samplesheet metadata/samplesheet.csv --genome data/assembly.fna -o results/sarek --more_args '--save_mapped --save_output_as_bam'"
    echo "  - To just print the help message for this script:"
    echo "      bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--samplesheet    <file>  Input samplesheet CSV file (see https://nf-co.re/sarek/3.2.3/usage)"
    echo "  --genome            <file>  Genome FASTA file"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --variant_tools     <str>   Comma-separated list of one or more tools for variant calling and annotation"
    echo "                              See https://nf-co.re/sarek/3.2.3/parameters#tools for options   [default: 'mpileup']"
    echo "  --download_wf       <str>   Download the specified workflow version (e.g. '--download_wf 3.2.0')"
    echo "                              The default current workflow version is $workflow_version"
    echo "  --workflow_dir      <dir>   Dir with the nf-core Sarek workflow                             [default: 'workflows/nfcore-sarek']"
    echo "                                - The workflow files will be downloaded here if needed."
    echo "  --more_args         <str>   Quoted string with more argument(s) for the workflow"
    echo "                              For a list of possibilities, see https://nf-co.re/sarek/3.2.3/parameters"
    echo
    echo "NEXTFLOW OPTIONS:"
    echo "  --restart                   Don't attempt to resume workflow run, but start over            [default: resume workflow]"
    echo "  --config            <file>  Additional config file                                          [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --profile           <str>   'Profile' to use from one of the config files                   [default: 'singularity']"
    echo "  --workdir           <dir>   Scratch (work) dir for the workflow                             [default: '/fs/scratch/PAS0471/\$USER/nfcore-sarek']"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v/--version                Print the version of this script and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
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
    # shellcheck source=/dev/null
    source "$function_script"
}

nextflow_env() {
    # Singularity container dir - any downloaded containers will be stored here;
    # if the required container is already there, it won't be re-downloaded
    export NXF_SINGULARITY_CACHEDIR="$container_dir"
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
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
samplesheet=
outdir=
genome=
more_args=
config_file=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --samplesheet )    shift && readonly samplesheet=$1 ;;
        -o | --outdir )         shift && readonly outdir=$1 ;;
        --genome )              shift && readonly genome=$1 ;;
        --variant_tools )       shift && readonly variant_tools=$1 ;;
        --workflow_dir )        shift && readonly workflow_dir=$1 ;;
        --download_wf )         shift && download_wf=true && workflow_version=$1 ;;
        --container_dir )       shift && readonly container_dir=$1 ;;
        --config )              shift && readonly config_file=$1 ;;
        --profile )             shift && readonly profile=$1 ;;
        --workdir )             shift && readonly work_dir=$1 ;;
        --restart )             resume=false ;;
        --more_args )           shift && readonly more_args=$1 ;;
        -v | --version )        script_version; exit 0 ;;
        -h | --help )           script_help; exit 0 ;;
        * )                     die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$samplesheet" ]] && die "No input samplesheet file specified, do so with -i/--samplesheet" "$all_args"
[[ -z "$genome" ]] && die "No input genome FASTA file specified, do so with --genome" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$samplesheet" ]] && die "Input samplesheet file $samplesheet does not exist"
[[ ! -f "$genome" ]] && die "Input genome FASTA file $genome does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict bash settings
set -euo pipefail

# Logging files and dirs
readonly LOG_DIR="$outdir"/logs
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
readonly TRACE_DIR="$outdir"/pipeline_info
mkdir -p "$LOG_DIR"

# Load software
load_env "$MODULE" "$CONDA" "$CONDA_YML"
nextflow_env

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Get the OSC config file
if [[ ! -f "$osc_config" ]]; then
    wget -q "$OSC_CONFIG_URL"
    osc_config=$(basename "$OSC_CONFIG_URL")
fi

# Build the config argument
config_arg="-c $osc_config"
if [[ "$config_file" != "" ]]; then
    config_arg="$config_arg -c ${config_file/,/ -c }"
fi

# Other
[[ "$resume" == false ]] && resume_arg=""
nf_file="$workflow_dir"/workflow/main.nf

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input samplesheet file:                   $samplesheet"
echo "Input genome FASTA file:                  $genome"
echo "Output dir:                               $outdir"
echo "Variant calling & annotation tools:       $variant_tools"
echo "Config file argument:                     $config_arg"
[[ -n $more_args ]] && echo "Additional arguments:                     $more_args"
log_time "Listing the input file(s):"
ls -lh "$samplesheet" "$genome"
log_time "Printing the contents of the samplesheet file:"
cat "$samplesheet"
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Make necessary dirs
mkdir -p "$work_dir" "$container_dir" "$TRACE_DIR"

# Download workflow, if needed
if [[ ! -d "$workflow_dir" || "$download_wf" == true ]]; then
    echo "# Downloading workflow to $workflow_dir"
    nf-core download sarek \
        --compress none \
        --container singularity \
        --outdir "$workflow_dir" \
        --revision "$workflow_version" \
        --force
    echo
elif [[ ! -f "$nf_file" ]]; then
    Die "Can't find workflow file $nf_file inside workflow dir $workflow_dir"
fi

# Run the tool
log_time "# Starting the Sarek workflow..."
runstats $TOOL_BINARY \
    "$nf_file" \
    --input "$samplesheet" \
    --outdir "$outdir" \
    --fasta "$genome" \
    --tools "$variant_tools" \
    $NONMODEL_OPTS \
    -work-dir "$work_dir" \
    -with-report "$TRACE_DIR"/report.html \
    -with-trace "$TRACE_DIR"/trace.txt \
    -with-timeline "$TRACE_DIR"/timeline.html \
    -with-dag "$TRACE_DIR"/dag.png \
    -ansi-log false \
    -profile "$profile" \
    $config_arg \
    $resume_arg \
    $more_args

# List the output
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*

# ==============================================================================
#                               WRAP UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
tool_version "$VERSION_COMMAND" | tee "$VERSION_FILE"
script_version "$SCRIPT_NAME" "$SCRIPT_AUTHOR" "$SCRIPT_VERSION" "$SCRIPT_URL" | tee -a "$VERSION_FILE" 
env | sort > "$ENV_FILE"
[[ "$IS_SLURM" = true ]] && resource_usage
log_time "Done with script $SCRIPT_NAME\n"
