#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=cutadapt
#SBATCH --output=slurm-cutadapt-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run Cutadapt to remove metabarcoding primers for a single pair of FASTQ files
The script will compute and use the reverse complements of all primers as well."
SCRIPT_VERSION="2025-02-28"
SCRIPT_AUTHOR="Jelmer Poelstra"
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=cutadapt
TOOL_NAME=Cutadapt
TOOL_DOCS=https://cutadapt.readthedocs.io/en/stable
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env_type=conda                     # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/cutadapt
container_url=
container_dir="$HOME/containers"

# Defaults - tool parameters
single_end=false
discard_untrimmed=true

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "\n                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage (always submit your scripts to SLURM with 'sbatch'):"
    echo "      sbatch $0 -i data/sample1_R1.fastq.gz -o results/cutadapt -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT"
    echo "  - Using a primer file instead:"
    echo "      sbatch $0 -i data/sample1_R1.fastq.gz -o results/cutadapt --primer_file metadata/primers.txt"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1             <file>  Input R1 FASTQ/FASTA file (name of R2 will be inferred in case of paired-end)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "There are two ways of specifying primers: (1) with '-f' and '-r' or (2) with '--primer_file' (use the latter if you have multiple primer pairs):"
    echo "  -f/--primer_f       <str>   Forward primer sequence (use in combination with -r)"
    echo "  -r/--primer_r       <str>   Reverse primer sequence (use in combination with -f)"
    echo "  --primer_file       <file>  File with primer sequences, one pair per line separated by a space (*alternative* to using -f and -r)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --single_end                Sequences are single-end                [default: $single_end]"
    echo "  --keep_untrimmed            Don't discard untrimmed sequences (i.e. those with no primers) [default: discard]"
    echo "  --more_opts         <str>   Quoted string with additional argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env_type          <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env_type]"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --container_path    <file>  Pre-existing Singularity container image file (.sif) to use"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "HARDCODED OPTIONS:"
    echo "  - The CutAdapt option '--pair-filter=any' is always used."
    echo 
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
    echo
}

# Function to source the script with Bash functions
source_function_script() {
    # Determine the location of this script, and based on that, the function script
    if [[ "$IS_SLURM" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script_name="$(basename "$FUNCTION_SCRIPT_URL")"
    function_script="$script_dir"/../dev/"$function_script_name"

    # Download the function script if needed, then source it
    if [[ -f "$function_script" ]]; then
        source "$function_script"
    elif [[ ! -f "$function_script_name" ]]; then
        echo "Can't find script with Bash functions ($function_script_name), downloading from GitHub..."
        wget -q "$FUNCTION_SCRIPT_URL" -O "$function_script_name"
        source "$function_script_name"
    else
        source "$function_script_name"
    fi
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
version_only=false  # When true, just print tool & script version info and exit
R1=
outdir=
discard_opt=
primer_f=
primer_r=
primer_file=
primer_opt=
more_opts=
threads=
container_path=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && R1=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        -f | --primer_f )   shift && primer_f=$1 ;;
        -r | --primer_r )   shift && primer_r=$1 ;;
        --primer_file )     shift && primer_file=$1 ;;
        --keep_untrimmed )  discard_untrimmed=false ;;
        --single_end )      single_end=true ;;
        --more_opts )       shift && more_opts=$1 ;;
        --env_type )        shift && env_type=$1 ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v | --version )         version_only=true ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Load software
load_env "$env_type" "$conda_path" "$container_path" "$container_dir" "$container_url"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$R1" ]] && die "No input file specified, do so with -i/--R1" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$R1" ]] && die "Input file $R1 does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
[[ "$discard_untrimmed" == true ]] && discard_opt="--discard-untrimmed"
indir=$(dirname "$R1")
R1_base=$(basename "$R1")
file_ext=$(echo "$R1_base" | sed -E 's/.*(.fastq|.fq|.fastq.gz|.fq.gz)$/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R1_basename=$(basename "$R1")
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")
[[ "$indir" == "$outdir" ]] && die "Input dir should not be the same as output dir ($indir)"

if [[ "$single_end" == false ]]; then
    log_time "Assuming paired-end sequences, inferring the R2 filename..."
    R2_suffix=${R1_suffix/1/2}
    R2="$indir"/${R1_base/$R1_suffix/$R2_suffix}
    R2_basename=$(basename "$R2")
    [[ ! -f "$R2" ]] && die "Input FASTQ file $R2 not found"
fi

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Working directory:                        $PWD"
echo "Input (R1) FASTQ file:                    $R1"
[[ "$single_end" == false ]] && echo "Input R2 FASTQ file:                      $R2"
echo "Output dir:                               $outdir"
echo "Discard untrimmed (-d):                   $discard_untrimmed"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$R1"
[[ "$single_end" == false ]] && ls -lh "$R2"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               DEFINE PRIMERS
# ==============================================================================
if [[ -z "$primer_file" ]]; then
    log_time "Using the forward and reverse primers provided as arguments to the script..."
    [[ -z "$primer_f" ]] && die "No forward primer (-f) provided"
    [[ -z "$primer_r" ]] && die "No reverse primer (-r) provided"

    primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
    primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

    if [[ "$single_end" == false ]]; then
        primer_opt="-a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"
    else
        primer_opt="-a $primer_f...$primer_r_rc"
    fi
else
    log_time "Using primer file $primer_file to read primers..."
    [[ ! -f "$primer_file" ]] && die "Primer file $primer_file not found"

    while read -r primer_f primer_r; do
        primer_f_rc=$(echo "$primer_f" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)
        primer_r_rc=$(echo "$primer_r" | tr ATCGYRKMBDHV TAGCRYMKVHDB | rev)

        if [[ "$single_end" == false ]]; then
            primer_opt="-a $primer_f...$primer_r_rc -A $primer_r...$primer_f_rc"
        else
            primer_opt="-a $primer_f...$primer_r_rc"
        fi
        # Remove leading whitespace:
        primer_opt=$(echo "$primer_opt" | sed -E 's/^ +//')
    done <"$primer_file"
fi

# Report
log_time "Forward primer (-f):              $primer_f"
log_time "Reverse primer (-r):              $primer_r"
log_time "Forward primer - rev. comp.:      $primer_f_rc"
log_time "Reverse primer - rev. comp.:      $primer_r_rc"
log_time "Primer option:                    $primer_opt"

# ==============================================================================
#                               RUN
# ==============================================================================
# Run Cutadapt
log_time "Running $TOOL_NAME..."
if [[ "$single_end" == false ]]; then
    runstats $TOOL_BINARY \
            $primer_opt \
            --output "$outdir"/"$R1_basename" \
            --paired-output "$outdir"/"$R2_basename" \
            --pair-filter=any \
            $discard_opt \
            --cores "$threads" \
            $more_opts \
            "$R1" "$R2"

    #? --pair-filter=any: Remove pair if one read is filtered (=Default)

else
    runstats $TOOL_BINARY \
            $primer_opt \
            --output "$outdir"/"$R1_basename" \
            $discard_opt \
            --cores "$threads" \
            $more_opts \
            "$R1"
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/"$sample_id"*
final_reporting "$LOG_DIR"
