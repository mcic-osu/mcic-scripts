#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=trimgalore
#SBATCH --output=slurm-trimgalore-%j.out

#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly DESCRIPTION="Run TrimGalore for one (single-end) or a pair of FASTQ files"
readonly SCRIPT_NAME=trimgalore.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly OSC_MODULE=miniconda3/4.12.0-py39
readonly CONDA_ENV=/fs/project/PAS0471/jelmer/conda/trimgalore
readonly TOOL_BINARY=trim_galore
readonly TOOL_NAME=TrimGalore
readonly TOOL_DOCS=https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

# Option defaults
quality=20                 # => 20 is also the TrimGalore default
length=20                  # => 20 is also the TrimGalore default
single_end=false           # => paired-end by default
nextseq=false              # => Assume non-2-color chemistry (used by NextSeq etc)
run_fastqc=true            # => Run FastQC after TrimGalore

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
    echo "      sbatch $0 -i data/fastq/S01_R1.fastq.gz -o results/trimgalore"
    echo "  - To run the script using a different OSC project than PAS0471:"
    echo "      sbatch -A PAS0001 $0 [...]"
    echo "  - To just print the help message for this script (-h) or for $TOOL_NAME (--help):"
    echo "      bash $0 -h"
    echo "      bash $0 --help"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1         <file>  Gzipped (R1) FASTQ input file (if paired-end, R2 file name will be inferred)"
    echo "  -o/--outdir     <dir>   Output dir"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  -q/--quality    <int>   Quality trimming threshold         [default: 20 (also the $TOOL_NAME default)]"
    echo "  -l/--length     <int>   Minimum read length                [default: 20 (also the $TOOL_NAME default)]"
    echo "  -s/--single_end         Input is single-end                [default: paired-end]"
    echo "  --nextseq               Sequences are from NextSeq or NovaSeq with 2-color chemistry"
    echo "                            This setting will help remove polyG tails"
    echo "  --no_fastqc             Don't run FastQC after trimming    [default: run FastQC after trimming]"
    echo "  --more_args             Additional arguments to pass to $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                      Print this help message and exit"
    echo "  --help                  Print the help for $TOOL_NAME and exit"
    echo "  -v                      Print the version of this script and exit"
    echo "  -v/--version            Print the version of $TOOL_NAME and exit"
    echo
    echo "OUTPUT:"
    echo "  Within the specified output dir (-o), results will be in 3 subdirectories:"
    echo "    - Trimmed FASTQ files will be in '<outdir>/trimmed/'"
    echo "    - FastQC output be in '<outdir>/fastqc/'"
    echo "    - Trimmomatic log files will be in '<outdir>/logs/'"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo
}

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
R1_in=""
outdir=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && readonly R1_in=$1 ;;
        -o | --outdir )     shift && readonly outdir=$1 ;;
        -q | --quality )    shift && readonly quality=$1 ;;
        -l | --length )     shift && readonly length=$1 ;;
        -s | --single_end ) readonly single_end=true ;;
        --nextseq )         readonly nextseq=true ;;
        --no_fastqc )       readonly run_fastqc=false ;;
        --more_args )       shift && readonly more_args=$1 ;;
        -v )                script_version; exit 0 ;;
        -h )                script_help; exit 0 ;;
        --version )         tool_version; exit 0 ;;
        --help )            tool_help; exit 0;;
        * )                 Die "Invalid option $1" "$all_args";;
    esac
    shift
done

# Check input
[[ -z "$R1_in" ]] && die "No input file specified, do so with -i/--R1" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ ! -f "$R1_in" ]] && die "Input file $R1_in does not exist"

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Source the Bash functions script
if [[ "$is_slurm" == true ]]; then
    SCRIPT_PATH=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
    SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
else
    SCRIPT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
fi
FUNCTION_SCRIPT="$SCRIPT_DIR"/../dev/bash_functions.sh
# shellcheck source=/dev/null
source "$FUNCTION_SCRIPT"

# Logging files and dirs
readonly log_dir="$outdir"/logs
readonly version_file="$log_dir"/version.txt
readonly conda_yml="$log_dir"/conda_env.yml
readonly env_file="$log_dir"/env.txt
mkdir -p "$log_dir"

# Load software and set nr of threads
load_tool_conda "$conda_yml"
set_threads

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Define outputs based on script parameters
outdir_trim="$outdir"/trimmed
outdir_fastqc="$outdir"/fastqc

# FastQC arguments
if [[ "$run_fastqc" == true ]]; then
    fastqc_arg1="--fastqc --fastqc_args"
    fastqc_arg2="-t $threads --outdir $outdir_fastqc"
else
    fastqc_arg1= && fastqc_arg2=
fi

# NextSeq arg
if [[ "$nextseq" == true ]]; then
    quality_arg="--nextseq $quality"
else
    quality_arg="--quality $quality"
fi

# Get file extension (.fastq.gz or .fq.gz)
extension=$(echo "$R1_in" | sed -E 's/.*(\.fa?s?t?q\.gz$)/\1/')

# Get R2 file, create input argument, define output files
if [ "$single_end" != "true" ]; then
    # Paired-end sequences
    R1_suffix=$(echo "$R1_in" | sed -E "s/.*(_R?1)_?[[:digit:]]*$extension/\1/")
    R2_suffix=${R1_suffix/1/2}
    R2_in=${R1_in/$R1_suffix/$R2_suffix}
    input_arg="--paired $R1_in $R2_in"
    
    [[ ! -f "$R2_in" ]] && Die "Input R2 FASTQ file $R2_in does not exist"
    [[ "$R1_in" == "$R2_in" ]] && Die "Input R1 and R2 FASTQ files are the same file"

    sample_id=$(basename "$R1_in" | sed -E "s/${R1_suffix}_?[[:digit:]]*${extension}//")
    R1_out="$outdir_trim"/"$sample_id"_R1.fastq.gz
    R2_out="$outdir_trim"/"$sample_id"_R2.fastq.gz
else
    # Single-end sequences
    input_arg="$R1_in"
    sample_id=$(basename "$R1_in" | sed "s/${extension}//")
    R1_out="$outdir_trim"/"$sample_id".fastq.gz
fi

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "R1 input file:                            $R1_in"
echo "Base output dir:                          $outdir"
echo "Sequence quality threshold:               $quality"
echo "Minimum sequence length:                  $length"
echo "Sequences are single-end:                 $single_end"
echo "Sequences are from NextSeq/NovaSeq:       $nextseq"
echo "Run FastQC:                               $run_fastqc"
echo "Number of threads/cores:                  $threads"
echo
[[ "$single_end" != "true" ]] && echo "R2 input file (inferred):                 $R2_in"
echo "Sample ID (inferred):                     $sample_id"
echo "Output dir - FastQC:                      $outdir_fastqc"
echo "R1 output file:                           $R1_out"
[[ "$single_end" != "true" ]] && echo "R2 output file:                           $R2_out"
[[ -n $more_args ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
log_time "Listing the input file(s):"
ls -lh "$R1_in"
[[ -n "$R2_in" ]] && ls -lh "$R2_in"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Make output dirs
mkdir -p "$outdir_trim" "$outdir_fastqc"

# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --output_dir "$outdir_trim" \
    $quality_arg \
    --length "$length" \
    --gzip \
    -j "$threads" \
    $more_args \
    $fastqc_arg1 \
    "$fastqc_arg2" \
    $input_arg

# Rename the output files
echo -e "\n# Moving output files..."
mv -v "$outdir_trim"/"$(basename "$R1_in")"_trimming_report.txt "$log_dir"

if [ "$single_end" != "true" ]; then
    mv -v "$outdir_trim"/"$sample_id"*_val_1.fq.gz "$R1_out"
    mv -v "$outdir_trim"/"$sample_id"*_val_2.fq.gz "$R2_out"
else
    R1_basename="$(basename "$R1_in" "$extension")"
    mv -v "$outdir_trim"/"$R1_basename"_trimmed.fq.gz "$R1_out"
fi

# List the output
log_time "Listing FASTQ output files:"
ls -lh "$R1_out"
[[ "$single_end" != "true" ]] && ls -lh "$R2_out"

# ==============================================================================
#                               WRAP UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
tool_version | tee "$version_file"
script_version | tee -a "$version_file"
env | sort > "$env_file"
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME\n"