#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=bactopia3_tools
#SBATCH --output=slurm-bactopia3_tools-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Bactopia v3 tools"
MODULE=miniconda3
CONDA=/fs/project/PAS0471/jelmer/conda/bactopia-dev
SCRIPT_VERSION="2023-07-20"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY=bactopia
TOOL_NAME="Bactopia Tools"
TOOL_DOCS=https://bactopia.github.io/
VERSION_COMMAND="bactopia --version"

# Constants - Nextflow and Bactopia generic settings
# Note: The samplesheet will be saved in $outdir/samplesheet.tsv 
OSC_CONFIG_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/nextflow/osc.config
OSC_CONFIG=mcic-scripts/nextflow/osc.config     # Will be downloaded if not present here
QUEUE_SIZE=100                                  # Nr of jobs to be submitted at once
MAX_CPUS=48
MAX_TIME=1440                                   # In hours
MAX_MEMORY=128                                  # In GB

CLEAN_WORKDIR_ARG="--cleanup_workdir"           # Always remove workdir after a successful run

# Defaults - Nextflow generics
work_dir_default=/fs/scratch/$SLURM_JOB_ACCOUNT/$USER/bactopia
container_dir=/fs/project/PAS0471/containers
profile=singularity
resume=true && resume_arg="-resume"

# Defaults - settings

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
    echo "      sbatch $0 -i results/bactopia -o results/bactopia/tools --tool pangenome --more_args '--use_roary'"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--bactopia_dir   <dir>   Dir with previous Bactopia results"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  --tool              <str>   Bactopia tool to run -- see https://bactopia.github.io/v2.2.0/bactopia-tools"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --restart                   Don't attempt to resume workflow run, but always start over     [default: resume any previous run]"
    echo 
    echo "UTILITY AND NEXTFLOW OPTIONS:"
    echo "  --work_dir          <dir>   Dir for initial workflow output files                           [default: /fs/scratch/$SLURM_JOB_ACCOUNT/$USER/bactopia]"
    echo "  --config            <file>  Additional Nextflow config file                                 [default: none - but mcic-scripts OSC config will be used]"
    echo "                                - Settings in this file will override default settings"
    echo "                                - Note that the mcic-scripts OSC config file will always be included, too"
    echo "                                  (https://github.com/mcic-osu/mcic-scripts/blob/main/nextflow/osc.config)"
    echo "  --container_dir     <dir>   Directory with/for stored container images                      [default: $container_dir]"
    echo "  --profile           <str>   Nextflow 'Profile' to use from one of the config files          [default: 'singularity']"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - $TOOL_DOCS"
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
    # Download the script if needed
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions.sh
    fi
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
indir=
outdir=
tool=
config_file=
work_dir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --bactopia_dir )   shift && indir=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --tool )                shift && tool=$1 ;;
        --container_dir )       shift && container_dir=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        --config )              shift && config_file=$1 ;;
        --profile )             shift && profile=$1 ;;
        --work_dir )            shift && work_dir=$1 ;;
        --restart )             resume=false && resume_arg= ;;
        -h | --help )           script_help; exit 0 ;;
        -v )                    script_version; exit 0 ;;
        --version )             load_env "$MODULE" "$CONDA"
                                tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                     die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$indir" ]] && die "No input dir specified, do so with -i/--bactopia_dir" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -d "$indir" ]] && die "Input dir $indir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict bash settings
set -euo pipefail

# Outputs - make paths absolute
[[ ! "$indir" =~ ^/ ]] && indir=$(realpath "$indir")
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"
[[ ! "$OSC_CONFIG" =~ ^/ ]] && OSC_CONFIG="$PWD"/"$OSC_CONFIG"

# Logging files and dirs
LOG_DIR="$outdir"/logs
VERSION_FILE="$LOG_DIR"/version.txt
CONDA_YML="$LOG_DIR"/conda_env.yml
ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software
load_env "$MODULE" "$CONDA" "$CONDA_YML"
nextflow_env
set_threads "$IS_SLURM"

# Get the OSC config file, then build the config argument
if [[ ! -f "$OSC_CONFIG" ]]; then
    wget -q "$OSC_CONFIG_URL"
    OSC_CONFIG=$(basename "$OSC_CONFIG_URL")
fi
[[ ! -f "$OSC_CONFIG" ]] && OSC_CONFIG="$outdir"/$(basename "$OSC_CONFIG_URL")
config_arg="-c $OSC_CONFIG"
[[ -n "$config_file" ]] && config_arg="$config_arg -c ${config_file/,/ -c }"

# Set Nextflow 'work' dir
if [[ -z "$work_dir" && "$IS_SLURM" == true ]]; then
    work_dir="${work_dir_default/pas/PAS}"  # 'pas' is in lowercase otherwise! 
else
    work_dir="$outdir"/work
fi

# Define outputs
outdir_base=$(dirname "$outdir")
run_name=$(basename "$outdir")

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo
echo "INPUT AND OUTPUT:"
echo "Input Bactopia dir:                       $indir"
echo "Output dir:                               $outdir"
[[ -n "$more_args" ]] && echo -e "\nSETTINGS:"
[[ -n "$more_args" ]] && echo "More args:                                $more_args"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run (if any):             $resume"
echo "Container dir:                            $container_dir"
echo "Scratch ('work') dir:                     $work_dir"
echo "Config 'profile':                         $profile"
echo "Config file argument:                     $config_arg"
[[ -n "$config_file" ]] && echo "Additional config file:                       $config_file"
echo "=========================================================================="
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Move into output dir
log_time "Changing into dir $outdir_base..."
cd "$outdir_base" || exit 1

# Run Bactopia
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --wf "$tool" \
    --bactopia "$indir" \
    --run_name "$run_name" \
    --outdir "$outdir" \
    $CLEAN_WORKDIR_ARG \
    --max_cpus $MAX_CPUS \
    --max_time $MAX_TIME \
    --max_memory $MAX_MEMORY \
    -qs $QUEUE_SIZE \
    --singularity_cache "$container_dir" \
    -profile "$profile" \
    -w "$work_dir" \
    $resume_arg \
    $config_arg \
    -ansi-log false \
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
