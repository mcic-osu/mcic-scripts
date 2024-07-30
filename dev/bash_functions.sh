#!/bin/bash

# Load software
load_env() {
    local module=$1
    local conda_env=$2
    local yml_file=${3-none}
    set +u

    # Load the OSC Conda module
    module load "$module"

    # Deactivate any active Conda environment
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi

    # Activate the focal environment
    log_time "Loading Conda environment $conda_env"
    source activate "$conda_env"

    # Store the Conda env in a YAML file
    [[ "$yml_file" != "none" ]] && conda env export --no-build > "$yml_file"

    set -u
}

# Print the tool's version
tool_version() {
    local version_command=${1-none}

    set +e
    
    if [[ "$version_command" == "none" ]]; then
        $TOOL_BINARY --version
    else
        eval $version_command
    fi
    
    set -e
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

# Print the script version
script_version() {
    local script_name=$1
    local script_author=$2
    local script_version=$3
    local script_URL=$4

    echo "Run using $script_name by $script_author, version $script_version ($script_URL)"
}

# Print SLURM job resource usage info
resource_usage() {
    echo
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
            log_time "Setting nr of threads based on SLURM_CPUS_PER_TASK to $threads"
        elif [[ -n "$SLURM_NTASKS" ]]; then
            readonly threads="$SLURM_NTASKS"
            log_time "Setting nr of threads based on SLURM_NTASKS to $threads"
        else 
            log_time "WARNING: Can't detect nr of threads, setting to 1"
            readonly threads=1
        fi
    else
        log_time "This is not a SLURM job, setting the nr of threads to 1"
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
    local version_command=$1
    local version_file=$2
    local env_file=$3
    local is_slurm=$4
    local script_name=$5
    local script_author=$6
    local script_version=$7
    local script_url=$8

    printf "\n======================================================================"
    log_time "Versions used:"
    tool_version "$version_command" | tee "$version_file"
    script_version "$script_name" "$script_author" "$script_version" "$script_url" | tee -a "$version_file" 
    env | sort > "$env_file"
    [[ "$is_slurm" = true ]] && resource_usage
    log_time "Done with script $script_name\n"
}
