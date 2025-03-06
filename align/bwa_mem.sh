#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=bwa
#SBATCH --output=slurm-bwa_mem-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
readonly DESCRIPTION="Map short reads to a reference using BWA MEM, and output a sorted BAM file and a 'samtools flagstat' QC stats file"
readonly MODULE=miniconda3
readonly CONDA=/fs/project/PAS0471/jelmer/conda/bwa-0.7.17
readonly SCRIPT_VERSION="2023-07-14"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY="bwa mem"
readonly TOOL_NAME=BWA
readonly TOOL_DOCS=https://github.com/lh3/bwa
readonly VERSION_COMMAND='bwa 2>&1 | grep Version'

# Option defaults
single_end=false            # Assume that reads are PE
use_secondary=true          # Use bwa mem's -M option: output additional shorter hits as secondary

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
    echo "  - Basic usage (always submit your scripts to SLURM with 'sbatch'):"
    echo "      sbatch $0 -i data/sampleA_R1.fastq.gz --index_dir results/bwa/index -o results/bwa"
    echo "  - To just print the help message for this script (-h) or for $TOOL_NAME (--help):"
    echo "      bash $0 -h"
    echo "      bash $0 --help"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  Input FASTQ file - for paired-end reads, provide only the R1"
    echo "                            (R2 file name will be inferred)"
    echo "  --index_dir     <dir>   Genome index dir"
    echo "  -o/--outdir     <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --no_secondary          Don't use BWA MEM's '-M' option to output shorter addtional hits as 'secondary' [default: use -M]"    
    echo "  --readgroup     <str>   Readgroup string to be added to the BAM file"
    echo "  --single_end            Reads are single-end - don't look for R2 file [default: PE reads]"
    echo "  --more_args     <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  --version               Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
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
    
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        git clone https://github.com/mcic-osu/mcic-scripts.git
        function_script=mcic-scripts/dev/bash_functions.sh
    fi
    # shellcheck source=/dev/null
    source "$function_script"
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
infile=
outdir=
index_dir=
R2=
more_args=
readgroup_string= && readgroup_arg=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && readonly infile=$1 ;;
        --index_dir )       shift && readonly index_dir=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        --single_end )      single_end=true ;;
        --no_secondary )    use_secondary=false ;;
        --readgroup )       shift && readonly readgroup_string=$1 ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v | --version )         load_env "$MODULE" "$CONDA"
                            tool_version "$VERSION_COMMAND" && exit 0 ;;
        * )                 die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check arguments
[[ -z "$infile" ]] && die "No input FASTQ file specified, do so with -i/--R1" "$all_args"
[[ -z "$index_dir" ]] && die "No input index dir specified, do so with --index_dir" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ ! -d "$index_dir" ]] && die "Input index dir $index_dir does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP II
# ==============================================================================
# Strict bash settings
set -euo pipefail

# Logging files and dirs
readonly LOG_DIR="$outdir"/logs
readonly VERSION_FILE="$LOG_DIR"/version.txt
readonly CONDA_YML="$LOG_DIR"/conda_env.yml
readonly ENV_FILE="$LOG_DIR"/env.txt
mkdir -p "$LOG_DIR"

# Load software and set nr of threads
load_env "$MODULE" "$CONDA" "$CONDA_YML"
set_threads "$IS_SLURM"

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Determine R2 file, output prefix, etc
if [[ -n "$infile" ]]; then
    R1_basename=$(basename "$infile" | sed -E 's/.fa?s?t?q.gz//')
    
    if [[ "$single_end" == false ]]; then
        R1_suffix=$(echo "$infile" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
        sampleID=${R1_basename/"$R1_suffix"/}

        # Determine name of R2 file
        if [[ -z "$R2" ]]; then
            R2_suffix=${R1_suffix/1/2}
            R2=${infile/$R1_suffix/$R2_suffix}
        fi
        
        [[ ! -f "$R2" ]] && die "R2 input file $R2 does not exist"
        [[ "$infile" == "$R2" ]] && die "Input file R1 is the same as R2 ($infile)"
    
    else
        sampleID="$R1_basename"
    fi
fi

# Reference index prefix
n_index=$(find "$index_dir" -name "*amb" | wc -l)
[[ "$n_index" == 0 ]] && die "No BWA index files found in index dir $index_dir"
[[ "$n_index" -gt 1 ]] && log_time "WARNING: More than one BWA indices found, using the first"
index_prefix=$(find "$index_dir" -name "*amb" | head -n 1 | sed 's/.amb//')

# Other
[[ -n "$readgroup_string" ]] && readgroup_arg="-R $readgroup_string"
if [[ "$use_secondary" == true ]]; then secondary_arg="-M"; else secondary_arg=; fi
bam=$outdir/"$sampleID".bam
flagstat_file="$outdir"/flagstat/"$sampleID".flagstat

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input (R1) FASTQ file:                    $infile"
[[ -n "$R2" ]] && echo "Input R2 FASTQ file:                      $R2"
echo "Reads are single end:                     $single_end"
echo "Shorter hits are output as secondary (-M option): $use_secondary"
echo "Index prefix:                             $index_prefix"
echo "Output BAM file:                          $bam"
echo "Output flagstat file:                     $flagstat_file"
[[ -n "$readgroup_string" ]] && echo "Readgroup string:                         $readgroup_string"
[[ -n $more_args ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ -n "$R2" ]] && ls -lh "$R2"
echo
ls -lh "${index_prefix}"*
[[ "$IS_SLURM" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Make output dirs
mkdir -p "$outdir"/flagstat

# Run BWA
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    -t "$threads" \
    $readgroup_arg \
    $secondary_arg \
    $more_args \
    "$index_prefix" \
    "$infile" "$R2" \
    | samtools sort \
        --threads "$threads" \
        -o "$bam" \
        -

# Get mapping stats
log_time "Getting mapping stats with Samtools flagstat:"
samtools flagstat "$bam" > "$flagstat_file"

# List the output
log_time "Listing the output BAM file:"
ls -lh "$bam"
log_time "Number of mapped reads:"
grep "primary mapped" "$flagstat_file"
[[ -n "$R2" ]] && grep "properly paired" "$flagstat_file" 

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

# ==============================================================================
#                               INFO
# ==============================================================================
#? Useful bwa flags:
#?  -t nr of threads
#?  -a alignments for single-end / unpaired reads are also output, as secondary alignments
#?  -M shorter split reads are output as secondary alignments, for Picard compatibility
#?  -R "@RG\tID:group1\tSM:$IND\tPL:illumina\tLB:lib1"
