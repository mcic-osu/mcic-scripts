#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=bactopia3
#SBATCH --output=slurm-bactopia3-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run the Bactopia v3 workflow"
MODULE=miniconda3
CONDA=/fs/project/PAS0471/jelmer/conda/bactopia-dev
SCRIPT_VERSION="2023-07-19"
SCRIPT_AUTHOR="Jelmer Poelstra"
SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
TOOL_BINARY=bactopia
TOOL_NAME=Bactopia
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
TRIM_ADAPTER_ARG="--trim"

# Defaults - Nextflow generics
work_dir_default=/fs/scratch/$SLURM_JOB_ACCOUNT/$USER/bactopia
container_dir=/fs/project/PAS0471/containers
profile=singularity
resume=true && resume_arg="-resume"

# Defaults - settings
use_bakta=true                                  # Differs from Bactopia default
assembler=spades                                # Differs from Bactopia default

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
    echo "      sbatch $0 -i data/fastq -o results/bactopia --species 'Salmonella enterica' --bakta_db /fs/ess/PAS0471/jelmer/refdata/bakta/db"
    echo "  - Use Prokka instead of Bakta for annotation:"
    echo "      sbatch $0 -i data/fastq -o results/bactopia --species 'Salmonella enterica' --use_prokka"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir          <dir>   Input dir with FASTQ files"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  --species           <str>   Focal species, e.g. 'Salmonella enterica' (make sure to quote!)"
    echo "  --bakta_db          <dir>   Dir with Bakta DB, required when using Bakta for annotation (use '--use_prokka' to use Prokka instead)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --assembler         <str>   Which assembler to use inside Shovill                           [default: spades]"
    echo "                              Options: spades, skesa, megahit, velvet"
    echo "  --use_prokka                Use Prokka instead of Bakta for annotation                      [default: Bakta]"
    echo "  --more_args         <str>   Quoted string with additional argument(s) for $TOOL_NAME"
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
    echo "                              Use 'none' to not load a profile at all"
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
species=
bakta_db=
bakta_arg=
config_file=
work_dir=
profile_arg=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )          shift && indir=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --species )             shift && species=$1 ;;
        --bakta_db )            shift && bakta_db=$1 ;;
        --use_prokka )          use_bakta=false ;;
        --assembler )           shift && assembler=$1 ;;
        --container_dir )       shift && container_dir=$1 ;;
        --more_args )           shift && more_args=$1 ;;
        --config )              shift && config_file=$1 ;;
        --profile )             shift && profile=$1 ;;
        --work_dir )            shift && work_dir=$1 ;;
        -restart )              resume=false && resume_arg= ;;
        -h | --help )           script_help; exit 0 ;;
        -v )                    script_version; exit 0 ;;
        --version )             load_env "$MODULE" "$CONDA"
                                tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                     die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$indir" ]] && die "No input dir specified, do so with -i/--indir" "$all_args"
[[ -z "$species" ]] && die "No species specified, do so with --species" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -d "$indir" ]] && die "Input dir $indir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict bash settings
set -euo pipefail

# Outputs - make paths absolute
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

# Profile
if [[ "$profile" != "none" ]]; then
    profile_arg="-profile $profile"
else
    [[ "$IS_SLURM" == true ]] && profile_arg="--queue slurm --cluster_opts '--account=PAS0471'"
fi

# Define outputs
#outdir_base=$(dirname "$outdir")
run_name=$(basename "$outdir") #TODO make option
outdir_main="$outdir"/main
outdir_cleanreads="$outdir"/cleanreads
fastqdir_cleanreads="$outdir_cleanreads"/fastq_out
samplesheet_raw="$outdir_cleanreads"/samplesheet.tsv
samplesheet_clean="$outdir_main"/samplesheet.tsv

# Bakta argument 
if [[ "$use_bakta" == true ]]; then
    [[ -z "$bakta_db" ]] && die "Using Bakta, but no Bakta DB dir specified, do so with --bakta_db" "$all_args"
    [[ ! -d "$bakta_db" ]] && die "Bakta DB dir $bakta_db does not exist"
    bakta_arg="--use_bakta --bakta_db $bakta_db"
fi

# Assembler arg
assembler_arg="--shovill_assembler $assembler"

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo
echo "INPUT AND OUTPUT:"
echo "Input FASTQ dir:                          $indir"
echo "Output dir:                               $outdir"
[[ "$use_bakta" == true ]] && echo "Bakta DB dir:                             $bakta_db"
echo "Nr of FASTQ files in the input dir:       $(ls "$indir"/*fastq.gz | wc -l)"
echo
echo "SETTINGS:"
echo "Species:                                  $species"
echo "Genome assembler (within Shovill):        $assembler"
echo "Using Bakta for annotation:               $use_bakta"
[[ -n "$more_args" ]] && echo "More args for Bactopia:                   $more_args"
echo
echo "NEXTFLOW-RELATED SETTINGS:"
echo "Resume previous run (if any):             $resume"
echo "Container dir:                            $container_dir"
echo "Scratch ('work') dir:                     $work_dir"
echo "Config 'profile' argument:                $profile_arg"
echo "Config file argument:                     $config_arg"
[[ -n "$config_file" ]] && echo "Additional config file:                       $config_file"
echo "=========================================================================="
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create output dirs
mkdir -p "$outdir_main" "$outdir_cleanreads" "$fastqdir_cleanreads"

# Move into output dir
#log_time "Changing into dir $outdir_base..."
#cd "$outdir_base" || exit 1

# Prepare the samplesheet
log_time "Preparing the samplesheet $samplesheet_raw with 'bactopia prepare'..."
runstats bactopia prepare \
    --path "$fastqdir_cleanreads" \
    --species "$species" \
    > "$samplesheet_raw"
log_time "Printing the contents of the sample sheet:"
cat -n "$samplesheet_raw"

# Run Bactopia - clean-yer-reads
log_time "Running the Bactopia 'clean-yer-reads' workflow..."
runstats $TOOL_BINARY \
    --wf cleanyerreads \
    --samples "$samplesheet_raw" \
    --outdir "$outdir_cleanreads" \
    --run_name "$run_name" \
    --species "$species" \
    $CLEAN_WORKDIR_ARG \
    --max_cpus $MAX_CPUS \
    --max_time $MAX_TIME \
    --max_memory $MAX_MEMORY \
    -qs $QUEUE_SIZE \
    --singularity_cache "$container_dir" \
    -w "$work_dir" \
    $profile_arg \
    $resume_arg \
    $config_arg \
    -ansi-log false

# Collect FASTQs
find "$outdir_cleanreads" -wholename "*/qc/*[12].fastq.gz" \
    -exec ln -sf {} "$fastqdir_cleanreads" \;

# Prepare the samplesheet
log_time "Preparing the samplesheet $samplesheet_clean with 'bactopia prepare'..."
runstats bactopia prepare \
    --path "$fastqdir_cleanreads" \
    --species "$species" \
    > "$samplesheet_clean"
log_time "Printing the contents of the sample sheet:"
cat -n "$samplesheet_clean"

# Run Bactopia - main workflow
log_time "Running the main Bactopia workflow..."
runstats $TOOL_BINARY \
    --wf bactopia \
    --samples "$samplesheet_clean" \
    --outdir "$outdir_main" \
    --run_name "$run_name" \
    --species "$species" \
    $assembler_arg \
    $bakta_arg \
    $CLEAN_WORKDIR_ARG \
    $TRIM_ADAPTER_ARG \
    --skip_qc \
    --max_cpus $MAX_CPUS \
    --max_time $MAX_TIME \
    --max_memory $MAX_MEMORY \
    -qs $QUEUE_SIZE \
    --singularity_cache "$container_dir" \
    -w "$work_dir" \
    $profile_arg \
    $resume_arg \
    $config_arg \
    -ansi-log false \
    $more_args

#? Other options
# --datasets_cache [string]  Directory where downloaded datasets should be stored. [default: <BACTOPIA_DIR>/data/datasets]
# --min_contig_len                    [integer] Minimum contig length <0=AUTO> [default: 500]
# --min_contig_cov                    [integer] Minimum contig coverage <0=AUTO> [default: 2]
# --organism                          [string]  Taxonomy group to run additional screens against (AMRFinder+) #TODO add this as an option

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
