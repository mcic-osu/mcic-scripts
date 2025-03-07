#!/bin/bash

# Constants
OSC_MODULE=miniconda3
REPO_URL=https://github.com/mcic-osu/mcic-scripts

# Dummy defaults
[[ -z "$env_type" ]] && env_type=conda
[[ -z "$container_path" ]] && container_path=
[[ -z "$SCRIPT_NAME" ]] && SCRIPT_NAME=script-name
[[ -z "$SCRIPT_VERSION" ]] && SCRIPT_VERSION=script-version
[[ -z "$SCRIPT_AUTHOR" ]] && SCRIPT_AUTHOR=script-author
[[ -z "$TOOL_NAME" ]] && TOOL_NAME=tool-name

# Variables that can/should be loaded in the script calling these functions
# conda_path        - Absolute path to a Conda environment dir
# container_url     - URL to a container
# container_path    - Absolute path to a container .sif file
# TOOL_BINARY       - The command that calls the focal program 
# SCRIPT_NAME       - Name of the shell script
# SCRIPT_VERSION    - Version of the shell script
# SCRIPT_AUTHOR     - Author of the shell script

# Load Conda or Singularity env
load_env() {
    if [[ "$env_type" == "conda" ]]; then
        load_conda
    elif [[ "$env_type" == "container" ]]; then
        load_container
    elif [[ "$env_type" == "none" ]]; then
        log_time "NOTE: not using a Conda environment OR a container, software expected to be in PATH"
    else
        die "Execution environment ('--env_type') should be 'conda', 'container', or 'none' but is currently $env_type"
    fi
}

# Load Conda env
load_conda() {
    set +u

    # Load the OSC Conda module
    module load "$OSC_MODULE"

    # Deactivate any active Conda environment
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do conda deactivate 2>/dev/null; done
    fi

    # Activate the focal environment
    log_time "Loading Conda environment $conda_path"
    conda activate "$conda_path"

    # No container prefix when using a Conda env
    CONTAINER_PREFIX=

    set -u
}

# Set up container
load_container() {
    dl_container=false
    
    # If no path to a container image file was provided,
    # then build the path based on the URL, and check if the file exists
    if [[ -z "$container_path" ]]; then
        url_basename=$(basename "$container_url")
        container_path="$container_dir"/${url_basename/:/_}.sif
        
        if [[ -f "$container_path" ]]; then
            log_time "No container path was provided, but the container in $container_url
            was found at $container_path and will be used"
        else
            dl_container=true
        fi
    fi

    # If needed, download the container image
    if [[ "$dl_container" == true ]]; then
        log_time "Downloading container from $container_url to $container_path"
        mkdir -p "$container_dir"
        singularity pull --force "$container_path" "$container_url"
    fi

    # Set the final 'prefix' to run the container
    CONTAINER_PREFIX="singularity exec $container_path"
    TOOL_BINARY="$CONTAINER_PREFIX $TOOL_BINARY"
    VERSION_COMMAND="$CONTAINER_PREFIX $VERSION_COMMAND"
    log_time "Using a container with base call: $CONTAINER_PREFIX"
}

# Print the tool's version
print_version() {
    local version_command=${1-none}
    set +e
    
    echo -e "\n# Version of this shell script:"
    echo "$SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION ($REPO_URL)"

    echo "# Version of $TOOL_NAME:"
    if [[ "$version_command" == "none" ]]; then
        $TOOL_BINARY --version
    else
        eval $CONTAINER_PREFIX $version_command
    fi
    
    set -e
}

# Print SLURM job resource usage info
resource_usage() {
    echo
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
}

# Print SLURM job requested resources
slurm_resources() {
    set +u
    log_time "SLURM job information:"
    echo "Account (project):                        $SLURM_JOB_ACCOUNT"
    echo "Job ID:                                   $SLURM_JOB_ID"
    echo "Job name:                                 $SLURM_JOB_NAME"
    echo "Memory (GB per node):                     $(( SLURM_MEM_PER_NODE / 1000 ))"
    echo "CPUs (on node):                           $SLURM_CPUS_ON_NODE"
    echo "Time limit (minutes):                     $(( SLURM_TIME_LIMIT / 60 ))"
    echo -e "==========================================================================\n"
    set -u
}

# Set the number of threads/CPUs
set_threads() {
    local is_slurm=$1
    set +u
    
    if [[ "$is_slurm" == true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            readonly threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            readonly threads="$SLURM_NTASKS"
        else 
            log_time "WARNING: This is a Slurm job but I can't detect the number of job threads, setting to 1"
            readonly threads=1
        fi
    else
        log_time "This is not a Slurm job, setting number of threads to 1"
        readonly threads=1
    fi
    
    export threads
    set -u
}

# Print command ran and its resource usage information for any process
runstats() {
    /usr/bin/time -f \
        "\n# Ran the command: \n%C
        \n# Run stats by /usr/bin/time:
        Time: %E   CPU: %P    Max mem: %M K    Exit status: %x \n" \
        "$@"
}

# Print log messages that include the time
log_time() {
    echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')]" ${1-""};
}

# Exit upon error with a message
die() {
    local error_message=${1}
    local error_args=${2-none}

    log_time "$0: ERROR: $error_message" >&2
    log_time "For help, run this script with the '-h' or '--help' option" >&2
    print_version "$VERSION_COMMAND"

    if [[ "$error_args" != "none" ]]; then
        log_time "Options passed to the script:" >&2
        echo "$error_args" >&2
    fi
    
    log_time "EXITING..." >&2
    exit 1
}

# Final reporting
final_reporting() {
    VERSION_FILE="$LOG_DIR"/versions.txt
    ENV_FILE="$LOG_DIR"/shell_env.txt
    CONDA_YML="$LOG_DIR"/conda_env.yml

    # Store the Conda env in a YAML file
    [[ "$env_type" == "conda" ]] && conda env export --no-build > "$CONDA_YML"

    printf "\n======================================================================"
    log_time "Versions used:"
    print_version "$VERSION_COMMAND" | tee "$VERSION_FILE" 
    env | sort > "$ENV_FILE"
    [[ "$IS_SLURM" == true ]] && resource_usage
    log_time "Done with script $SCRIPT_NAME\n"
}
