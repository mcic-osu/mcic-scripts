#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=nfc_ampliseq
#SBATCH --output=slurm-nfc_ampliseq-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run the Nextflow-core metabarcoding pipeline from https://nf-co.re/ampliseq"
SCRIPT_VERSION="2025-02-19"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY="nextflow run"
export TOOL_NAME="nextflow"
VERSION_COMMAND="nextflow -v"

# Constants - parameters
WORKFLOW_NAME=nf-core/ampliseq                              # The name of the nf-core workflow
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config

# Parameter defaults - infrastructure
conda_path=/fs/ess/PAS0471/jelmer/conda/nextflow
osc_account=PAS0471                                         # If the script is submitted with another project, this will be updated (line below)
[[ -n $SLURM_JOB_ACCOUNT ]] && osc_account=$(echo "$SLURM_JOB_ACCOUNT" | tr "[:lower:]" "[:upper:]")

# Parameter defaults - workflow
workflow_version=2.12.0                                     # The version of the nf-core workflow
work_dir=/fs/scratch/"$osc_account"/$USER/nfc-ampliseq      # 'work dir' for initial outputs (selected, final outputs go to the outdir)
profile="singularity"
resume=true && resume_opt="-resume"
container_dir="$work_dir"/containers                        # The workflow will download containers to this dir

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo "  - Different from the Nextflow default, this script will try to 'resume' (rather than restart) a previous incomplete run by default"
    echo "  - This workflow can be used for both 16S and ITS data: default is 16S; change the settings in nfcore_ampliseq_params.yml for ITS"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  sbatch $0 -o results/ampliseq -p config/nfcore_ampliseq_params.yml"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outdir         <dir>   Dir for pipeline output files (will be created if needed)"
    echo "  -p/--params         <file>  YAML file with workflow parameters"
    echo "                              Use 'mcic-scripts/metabar/nfcore_ampliseq_params.yml' as template"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --workflow_version  <str>   Nf-core ampliseq workflow version to use                    [default: $workflow_version]"
    echo "  --restart                   Don't attempt to resume workflow run, but start over        [default: resume workflow]"
    echo
    echo "ADVANCED NEXTFLOW-RELATED OPTIONS:"
    echo "  --work_dir           <dir>  Scratch (work) dir for the workflow                         [default: $work_dir]"
    echo "                                - This is where the workflow results will be stored before final results are copied to the output dir."
    echo "  --container_dir     <dir>   Directory with container images                             [default: $container_dir]"
    echo "                                - Required images will be downloaded here when not already present here" 
    echo "  --config            <file>  Additional config file                                      [default: none]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --profile            <str>  'Profile' to use from one of the config files               [default: $profile]"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of Nextflow and exit"
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
version_only=false                    # When true, just print tool & script version info and exit
outdir=
params_file=
config_file=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )                 shift && outdir=$1 ;;
        -p | --params )                 shift && params_file=$1 ;;
        --workflow_version )            shift && workflow_version=$1 ;;
        --container_dir )               shift && container_dir=$1 ;;
        --config | -config )            shift && config_file=$1 ;;
        --profile | -profile )          shift && profile=$1 ;;
        --work_dir | -work-dir )        shift && work_dir=$1 ;;
        --restart | -restart )          resume=false && resume_opt= ;;
        -h | --help )                   script_help; exit 0 ;;
        -v )                            script_version; exit 0 ;;
        --version )                     version_only=true ;;
        * )                             die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# Check arguments
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ -z "$params_file" ]] && die "No parameter YAML file specified, do so with -p/--params" "$all_opts"
[[ ! -f "$params_file"  ]] && die "Input parameter YAML file $params_file does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Load software and set nr of threads
load_env "$conda_path"
nextflow_setup
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Build the config argument
OSC_CONFIG="$outdir"/$(basename "$OSC_CONFIG_URL")
config_opt="-c $OSC_CONFIG"
[[ -n "$config_file" ]] && config_opt="$config_opt -c ${config_file/,/ -c }"

# Other dirs
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
trace_dir="$outdir"/pipeline_info

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:                 $all_opts"
echo
echo "INPUT AND OUTPUT:"
echo "Parameter YAML file:                          $params_file"
echo "Output dir:                                   $outdir"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run (if any):                 $resume"
echo "Container dir:                                $container_dir"
echo "Scratch (work) dir:                           $work_dir"
echo "Config 'profile':                             $profile"
echo "Config file argument:                         $config_opt"
[[ -n "$config_file" ]] && echo "Additional config file:                       $config_file"
echo "=========================================================================="
set_threads "$IS_SLURM"
[[ "$IS_SLURM" = true ]] && slurm_resources
echo "=========================================================================="
log_time "Printing the contents of the parameter file:"
cat -n "$params_file"
if [[ -n "$config_file" ]]; then
    log_time "Printing the contents of the additional config file:"
    cat -n "$config_file"
fi
echo "=========================================================================="

# ==============================================================================
#                               RUN
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

# Run the tool
log_time "Starting the workflow.."
runstats $TOOL_BINARY $WORKFLOW_NAME \
    -r "$workflow_version" \
    -params-file "$params_file" \
    --outdir "$outdir" \
    -work-dir "$work_dir" \
    -with-report "$trace_dir"/report.html \
    -with-trace "$trace_dir"/trace.txt \
    -with-timeline "$trace_dir"/timeline.html \
    -with-dag "$trace_dir"/dag.png \
    -ansi-log false \
    -profile "$profile" \
    $config_opt \
    $resume_opt

# List the output, report version, etc
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
