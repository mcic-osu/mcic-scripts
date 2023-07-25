#!/bin/bash

# Constants
OSC_MODULE=miniconda3

# Load Conda or Singularity env
load_env() {
    if [[ "$env" == "conda" ]]; then
        load_conda
    elif [[ "$env" == "container" ]]; then
        load_container
    else
        die "Execution environment ('--env') should be 'conda' or 'container' but is currently $env"
    fi
}

# Load Conda env
load_conda() {
    set +u

    # Load the OSC Conda module
    module load "$OSC_MODULE"

    # Deactivate any active Conda environment
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi

    # Activate the focal environment
    log_time "Loading Conda environment $conda_path"
    source activate "$conda_path"

    # No container prefix when using a Conda env
    CONTAINER_PREFIX=

    set -u
}

# Set up container
load_container() {
    if [[ "$dl_container" == true ]]; then
        log_time "Downloading container from $container_url to $container_dir"
        mkdir -p "$container_dir"
        singularity pull --dir "$container_dir" "$container_url"
        container_path="$container_dir"/$(basename "$container_url")
    fi

    CONTAINER_PREFIX="singularity exec $container_path"
    log_time "Using a container with base call: $CONTAINER_PREFIX"
}

# Print the tool's version
tool_version() {
    local version_command=${1-none}
    set +e
    
    if [[ "$version_command" == "none" ]]; then
        $CONTAINER_PREFIX $TOOL_BINARY --version
    else
        eval $CONTAINER_PREFIX $version_command
    fi
    
    set -e
}

# Print the script version
script_version() {
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
    echo "Account (project):              $SLURM_JOB_ACCOUNT"
    echo "Job ID:                         $SLURM_JOB_ID"
    echo "Job name:                       $SLURM_JOB_NAME"
    echo "Memory (GB per node):           $(( SLURM_MEM_PER_NODE / 1000 ))"
    echo "CPUs (per task):                $SLURM_CPUS_PER_TASK"
    echo "Time limit (minutes):           $(( SLURM_TIME_LIMIT / 60 ))"
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
        log_time "Arguments passed to the script:" >&2
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
    VERSION_FILE="$LOG_DIR"/version.txt
    ENV_FILE="$LOG_DIR"/env.txt
    CONDA_YML="$LOG_DIR"/conda_env.yml

    # Store the Conda env in a YAML file
    [[ "$env" == "conda" ]] && conda env export --no-build > "$CONDA_YML"

    printf "\n======================================================================"
    log_time "Versions used:"
    tool_version "$VERSION_COMMAND" | tee "$VERSION_FILE"
    script_version | tee -a "$VERSION_FILE" 
    env | sort > "$ENV_FILE"
    [[ "$IS_SLURM" = true ]] && resource_usage
    log_time "Done with script $SCRIPT_NAME\n"
}
