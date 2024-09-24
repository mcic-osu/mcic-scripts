#!/bin/bash

# Constants
OSC_MODULE=miniconda3

# Defaults
[[ -z "$env" ]] && env=conda
[[ -z "$REPO_URL" ]] && REPO_URL=repo-URL
[[ -z "$SCRIPT_NAME" ]] && SCRIPT_NAME=script-name
[[ -z "$SCRIPT_VERSION" ]] && SCRIPT_VERSION=script-version
[[ -z "$SCRIPT_AUTHOR" ]] && SCRIPT_AUTHOR=script-author
[[ -z "$TOOL_NAME" ]] && TOOL_NAME=tool-name

# Variables that should be loaded in the script calling these functions
# conda_path        Absolute path to a Conda environment dir
# container_url
# container_path
# TOOL_BINARY
# REPO_URL
# SCRIPT_NAME
# SCRIPT_VERSION
# SCRIPT_AUTHOR

# Load Conda or Singularity env
load_env() {
    if [[ "$env" == "conda" ]]; then
        load_conda
    elif [[ "$env" == "container" ]]; then
        load_container
    elif [[ "$env" == "none" ]]; then
        log_time "NOTE: not using a Conda environment OR a container, software expected to be in PATH"
    else
        die "Execution environment ('--env') should be 'conda', 'container', or 'none' but is currently $env"
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
    log_time "Using a container with base call: $CONTAINER_PREFIX"
}

# Print the tool's version
tool_version() {
    local version_command=${1-none}
    set +e
    
    echo "# Version of $TOOL_NAME:"
    if [[ "$version_command" == "none" ]]; then
        $CONTAINER_PREFIX $TOOL_BINARY --version
    else
        eval $CONTAINER_PREFIX $version_command
    fi
    
    set -e
}

# Print the script version
script_version() {
    echo "# Version of this shell script:"
    echo "$SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION ($REPO_URL)"
}

# Print the tool's help
tool_help() {
    local help_command=${1-none}

    if [[ "$help_command" == "none" ]]; then
        $TOOL_BINARY --help
    else
        eval $help_command
    fi
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
            #log_time "Setting nr of threads (based on SLURM_CPUS_PER_TASK) to $threads"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            readonly threads="$SLURM_NTASKS"
            #log_time "Setting nr of threads (based on SLURM_NTASKS) to $threads"
        else 
            log_time "WARNING: Can't detect nr of Slurm job threads, setting to 1"
            readonly threads=1
        fi
    else
        log_time "This is not a Slurm job, setting nr of threads to 1"
        readonly threads=1
    fi
    
    export threads
    set -u
}

# Resource usage information for any process
runstats() {
    /usr/bin/time -f \
        "\n# Ran the command: \n%C
        \n# Run stats by /usr/bin/time:
        Time: %E   CPU: %P    Max mem: %M K    Exit status: %x \n" \
        "$@"
}

# Exit upon error with a message
die() {
    local error_message=${1}
    local error_args=${2-none}

    log_time "$0: ERROR: $error_message" >&2
    log_time "For help, run this script with the '-h' option" >&2
    
    if [[ "$error_args" != "none" ]]; then
        log_time "Options passed to the script:" >&2
        echo "$error_args" >&2
    fi
    
    log_time "EXITING..." >&2
    exit 1
}

# Log messages that include the time
log_time() {
    echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')]" ${1-""};
}

# Final reporting
final_reporting() {
    VERSION_FILE="$LOG_DIR"/versions.txt
    ENV_FILE="$LOG_DIR"/shell_env.txt
    CONDA_YML="$LOG_DIR"/conda_env.yml

    # Store the Conda env in a YAML file
    [[ "$env" == "conda" ]] && conda env export --no-build > "$CONDA_YML"

    printf "\n======================================================================"
    log_time "Versions used:"
    tool_version "$VERSION_COMMAND"
    script_version
    
    script_version > "$VERSION_FILE"
    tool_version "$VERSION_COMMAND" &>> "$VERSION_FILE"
    env | sort > "$ENV_FILE"
    
    [[ "$IS_SLURM" == true ]] && resource_usage
    
    log_time "Done with script $SCRIPT_NAME\n"
}
