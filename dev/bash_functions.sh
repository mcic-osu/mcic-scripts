#!/bin/bash

# Load software
load_tool_conda() {
    local yml_file=${2-none}
    set +u

    # Load the OSC Conda module
    module load "$OSC_MODULE" 
    # Deactivate any active Conda environment
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi

    # Activate the focal environment
    source activate "$CONDA_ENV"

    # Store the Conda env in a YAML file
    [[ "$yml_file" != "none" ]] && conda env export --no-build > "$yml_file"

    set -u
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

# Print the script version
script_version() {
    echo "Run using $SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION ($SCRIPT_URL)"
}

# Print the tool's version
tool_version() {
    local version_command=${2-none}

    set +e
    load_tool_conda
    
    if [[ "$version_command" != none ]]; then
        $TOOL_BINARY --version
    else
        $version_command
    fi
    
    set -e
}

# Print the tool's help
tool_help() {
    local help_command=${2-none}

    load_tool_conda
    
    if [[ "$help_command" != none ]]; then
        $TOOL_BINARY --help
    else
        $help_command
    fi
}

# Print SLURM job resource usage info
resource_usage() {
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime |
        grep -Ev "ba|ex"
}

# Print SLURM job requested resources
slurm_resources() {
    set +u
    log_time "SLURM job information:"
    echo "Account (project):              $SLURM_JOB_ACCOUNT"
    echo "Job ID:                         $SLURM_JOB_ID"
    echo "Job name:                       $SLURM_JOB_NAME"
    echo "Memory (MB per node):           $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):                $SLURM_CPUS_PER_TASK"
    echo "Time limit:                     $SLURM_TIME_LIMIT"
    echo -e "=================================================================\n"
    set -u
}

# Set the number of threads/CPUs
set_threads() {
    set +u
    
    if [[ "$is_slurm" == true ]]; then
        if [[ -n "$SLURM_CPUS_PER_TASK" ]]; then
            readonly threads="$SLURM_CPUS_PER_TASK"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            readonly threads="$SLURM_NTASKS"
        else 
            log_time "WARNING: Can't detect nr of threads, setting to 1"
            readonly threads=1
        fi
    else
        readonly threads=1
    fi
    
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
