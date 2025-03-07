#!/usr/bin/env bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=hisat2
#SBATCH --output=slurm-hisat2-%j.out

# Run HISAT2 for DNA or RNA alignment to a reference genome

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=hisat2.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/hisat2
readonly TOOL_BINARY=hisat2
readonly TOOL_NAME=HISAT2
readonly TOOL_DOCS=http://daehwankimlab.github.io/hisat2/manual/
readonly TOOL_PAPER=https://www.nature.com/articles/s41587-019-0201-4

# Option defaults
single_end=false
#> The default strand direction is reverse

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  Run HISAT2 for DNA or RNA alignment to a reference genome"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  sbatch $0 --reads data/fastq/sampleA_R1.fastq.gz --index results/hisat2/index -o results/hisat2"
    echo "  sbatch $0 --reads data/fastq/sampleA.fastq.gz --index results/hisat2/index -o results/hisat2 --single_end"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --reads         <file>  (Optionally gzipped) FASTQ file"
    echo "                            For paired-end reads, provide only R1, the R2 filename will be inferred"
    echo "  --index_dir  <file>     Directory with the genome index"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --strandedness  <str>   Read strandedness: RF, FR, R, or F"
    echo "                            [default: RF for paired-end, R for single-end reads]"
    echo "  --more_args     <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for $TOOL_NAME and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
    echo
}

# Load software
load_tool_conda() {
    set +u
    module load "$MODULE" # Load the OSC Conda module
    # Deactivate any active Conda environments:
    if [[ -n "$CONDA_SHLVL" ]]; then
        for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    fi
    source activate "$CONDA_ENV" # Activate the focal environment
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
log_time() { echo -e "\n[$(date +'%Y-%m-%d %H:%M:%S')]" ${1-""}; }

# Print the script version
script_version() {
    echo "Run using $SCRIPT_NAME by $SCRIPT_AUTHOR, version $SCRIPT_VERSION ($SCRIPT_URL)"
}

# Print the tool's version
print_version() {
    set +e
    load_tool_conda
    "$TOOL_BINARY" --version
    set -e
}

# Print the tool's help
tool_help() {
    load_tool_conda
    "$TOOL_BINARY" --help
}

# Print SLURM job resource usage info
resource_usage() {
    sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%60,Elapsed,CPUTime | grep -Ev "ba|ex"
}

# Print SLURM job requested resources
slurm_resources() {
    set +u
    log_time "SLURM job information:"
    echo "Account (project):                        $SLURM_JOB_ACCOUNT"
    echo "Job ID:                                   $SLURM_JOB_ID"
    echo "Job name:                                 $SLURM_JOB_NAME"
    echo "Memory (MB per node):                     $SLURM_MEM_PER_NODE"
    echo "CPUs (per task):                          $SLURM_CPUS_PER_TASK"
    echo "Time limit:                               $SLURM_TIMELIMIT"
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

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
reads=
outdir=
strandedness=
more_args=
index_dir=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --reads )           shift && readonly reads=$1 ;;
        --index_dir )       shift && readonly index_dir=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        --single_end )      shift && readonly single_end=true ;;
        --strandedness )    shift && readonly strandedness=$1 ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -h )                script_help; exit 0 ;;
        -v | --version )         print_version; exit 0 ;;
        --help )            tool_help; exit 0;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$reads" ]] && die "No input FASTQ file specified, do so with --reads" "$all_args"
[[ -z "$index_dir" ]] && die "No genome index dir specified, do so with --index_dir" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$reads" ]] && die "Input file $reads does not exist"
[[ ! -d "$index_dir" ]] && die "Genome index dir $index_dir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
load_tool_conda
set_threads

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Define outputs based on script parameters
readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs

# Determine R2 file, output prefix, etc
R1_basename=$(basename "$reads" | sed -E 's/.fa?s?t?q.gz//')
    
if [ "$single_end" = false ]; then
    R1_suffix=$(echo "$reads" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
    sampleID=${R1_basename/"$R1_suffix"/}

    R2_suffix=${R1_suffix/1/2}
    reads_R2=${reads/$R1_suffix/$R2_suffix}
    [[ ! -f "$reads_R2" ]] && die "R2 input file $reads_R2 does not exist"
    [[ "$reads" == "$reads_R2" ]] && die "Input file R1 is the same as R2: $reads"

    reads_arg="-1 $reads -2 $reads_R2"
else
    sampleID="$R1_basename"
    reads_arg="-U $reads"
fi

# Strandedness argument
if [[ "$single_end" == false ]]; then
    [[ -z "$strandedness" ]] && strand_arg="--rna-strandness RF"
    [[ -n "$strandedness" ]] && strand_arg="--rna-strandness $strandedness"
else 
    [[ -z "$strandedness" ]] && strand_arg="--rna-strandness R"
    [[ -n "$strandedness" ]] && strand_arg="--rna-strandness $strandedness"
fi

# Define other outputs
out_prefix="$outdir/$sampleID"
index_prefix=$(ls -1 "$index_dir"/*ht2 | head -n 1 | sed -E 's/\.1.ht2$//')

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Output dir:                               $outdir"
echo "Input (R1) FASTQ:                         $reads"
[[ "$single_end" == false ]] && echo "Input R2 FASTQ:                           $reads_R2"
echo "Genome index dir:                         $index_dir"
echo "Strandedness argument:                    $strand_arg"
[[ -n $more_args ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
echo "Number of threads/cores:                  $threads"
echo
log_time "Listing the input file(s):"
ls -lh "$reads" 
[[ "$single_end" == false ]] && ls -lh "$reads_R2"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir"

# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    -x "$index_prefix" \
    $reads_arg \
    $strand_arg \
    --summary-file "$out_prefix".hisat2.summary.log \
    --threads "$threads" \
    $more_args \
    | samtools view -bS -F 4 -F 256 - > "$out_prefix".bam

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
print_version | tee "$version_file"
script_version | tee -a "$version_file" 
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME"
echo
