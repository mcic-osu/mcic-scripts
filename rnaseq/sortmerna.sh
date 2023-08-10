#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sortmerna
#SBATCH --output=slurm-sortmerna-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Run SortMeRNA to separate paired-end RNAseq reads into rRNA-derived and other reads
Separate pairs of FASTQ files with reads classified as rRNA (mapped) and reads not classified
as rRNA (unmapped)"
SCRIPT_VERSION="2023-08-10"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY=sortmerna
TOOL_NAME=SortMeRNA
TOOL_DOCS=https://github.com/biocore/sortmerna
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/sortmerna-env
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true

# Constants - tool parameters
# Paths to rRNA reference files within the SortMeRNA repo dir
R18S_PATH=data/rRNA_databases/silva-euk-18s-id95.fasta  
R28S_PATH=data/rRNA_databases/silva-euk-28s-id98.fasta

# Defaults - tool parameters
deinterleave=true

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo "                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Basic usage:"
    echo "      sbatch $0 data/S1_R1.fastq.gz -o results/sortmerna"
    echo "  - Loop:"
    echo "      for R1 in data/fastq/*R1.fastq.gz; do"
    echo "          sbatch $0 -i \$R1 -o results/sortmerna"
    echo "      done"
    echo
    echo "OUTPUT:"
    echo "  - Directory '<outdir>/by_sample/<sample_id> will contain e.g. SortMeRNA DB files and log files"
    echo "  - Directories '<outdir>/mapped' and '<outdir>/unmapped' will contain FASTQ files (for all samples)" 
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--R1             <file>  Input R1 FASTQ file (name of R2 will be inferred)"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "                                When running this script for multiple samples (looping over samples),"
    echo "                                you can specify a single output dir for all of them (see example above)." 
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --repo              <dir>   Directory with SortMeRNA repo (with rRNA db) [default: download repo]"
    echo "  --as_is                     Don't 'de-interleave' output FASTQ file  [default: de-interleave]"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "                                (NOTE: If no default '--container_url' is listed below,"
    echo "                                 you'll have to provide one in order to run the script with a container.)"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or '--dl_container' is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: $dl_container]"
    echo "  --no_strict                 Don't use strict Bash settings ('set -euo pipefail') -- can be useful for troubleshooting"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
    echo
    echo "TOOL DOCUMENTATION: $TOOL_DOCS"
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
    function_script=$(realpath "$script_dir"/../dev/"$(basename "$FUNCTION_SCRIPT_URL")")
    # Download the function script if needed, then source it
    if [[ ! -f "$function_script" ]]; then
        echo "Can't find script with Bash functions ($function_script), downloading from GitHub..."
        function_script=$(basename "$FUNCTION_SCRIPT_URL")
        wget "$FUNCTION_SCRIPT_URL" -O "$function_script"
    fi
    source "$function_script"
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "$SLURM_JOB_ID" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
R1=
outdir=
opts=
version_only=false
threads=
repo_dir=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --R1 )         shift && R1=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --as_is )           deinterleave=false ;;
        --repo )            repo_dir=false ;;
        --opts )            shift && opts=$1 ;;
        --env )             shift && env=$1 ;;
        --no_strict )       strict_bash=false ;;
        --dl_container )    dl_container=true ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 && dl_container=true ;;
        -h | --help )       script_help; exit 0 ;;
        -v )                script_version; exit 0 ;;
        --version )         version_only=true ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
[[ "$strict_bash" == true ]] && set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$R1" ]] && die "No input file specified, do so with -i/--R1" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$R1" ]] && die "Input file $R1 does not exist"

# Define outputs based on script parameters
# Infer the name of the R2 file
file_ext=$(basename "$R1" | sed -E 's/.*(.fastq.gz|.fq.gz)/\1/')
R1_suffix=$(basename "$R1" "$file_ext" | sed -E "s/.*(_R?1)_?[[:digit:]]*/\1/")
R2_suffix=${R1_suffix/1/2}
R2=${R1/$R1_suffix/$R2_suffix}
sample_id=$(basename "$R1" "$file_ext" | sed -E "s/${R1_suffix}_?[[:digit:]]*//")
[[ ! -f "$R2" ]] && die "Input file $R2 does not exist"
[[ "$R1" == "$R2" ]] && die "Input file $R1 and $R2 are the same ($R1)"

# Define output files
outdir_sample="$outdir"/by_sample/"$sample_id"
LOG_DIR="$outdir_sample"/logs
out_mapped_raw="$outdir_sample"/mapped_raw/"$sample_id"
out_unmapped_raw="$outdir_sample"/unmapped_raw/"$sample_id"

R1_mapped="$outdir"/mapped/"$sample_id""$R1_suffix".fastq.gz
R2_mapped="$outdir"/mapped/"$sample_id""$R2_suffix".fastq.gz
R1_unmapped="$outdir"/unmapped/"$sample_id""$R1_suffix".fastq.gz
R2_unmapped="$outdir"/unmapped/"$sample_id""$R2_suffix".fastq.gz

# Reference FASTA files (to be downloaded)
[[ -z $repo_dir ]] && repo_dir="$outdir_sample"/sortmerna_repo
ref_18s="$repo_dir"/"$R18S_PATH"
ref_28s="$repo_dir"/"$R28S_PATH"

# Make output dirs
mkdir -p "$LOG_DIR" "$repo_dir" \
    "$outdir"/mapped "$outdir"/unmapped \
    "$outdir_sample"/mapped_raw "$outdir_sample"/unmapped_raw 

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "R1 input FASTQ:                           $R1"
echo "R2 input FASTQ:                           $R2"
echo "Output dir:                               $outdir"
echo "SortMeRNA repo dir:                       $repo_dir"
echo "Deinterleave FASTQ files:                 $deinterleave"
echo "18S reference file:                       $ref_18s"
echo "28S reference file:                       $ref_28s"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:        $opts"
log_time "Listing the input file(s):"
ls -lh "$R1" "$R2"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Clone sortmerna repo to get db FASTA files
if [[ ! -f "$ref_18s" && ! -f "$ref_28s" ]]; then
    log_time "Cloning SortMeRNA repo..."
    n_seconds=$(( RANDOM % 50 + 1 ))
    sleep "$n_seconds"s # Sleep for a while so git doesn't error when running this multiple times in parallel
    [[ ! -f "$ref_18s" && ! -f "$ref_28s" ]] && git clone https://github.com/biocore/sortmerna "$repo_dir"
fi
# Check that db files are there
[[ ! -f "$ref_18s" ]] && die "18s reference FASTA file $ref_18s not found"
[[ ! -f "$ref_28s" ]] && die "28s reference FASTA file $ref_28s not found"

# Run SortMeRNA
log_time "Running $TOOL_NAME..."
runstats $CONTAINER_PREFIX $TOOL_BINARY \
    --ref "$ref_18s" \
    --ref "$ref_28s" \
    --reads "$R1" \
    --reads "$R2" \
    --fastx \
    --aligned "$out_mapped_raw" \
    --other "$out_unmapped_raw" \
    --workdir "$outdir_sample" \
    --paired_in \
    --threads "$threads" \
    $opts

#?--paired_in Flags the paired-end reads as Aligned, when either of them is Aligned.

# De-interleave the output
if [[ "$deinterleave" = true ]]; then
    log_time "Deinterleaving mapped reads..."
    reformat.sh \
        in="$out_mapped_raw".fq.gz \
        out1="$R1_mapped" out2="$R2_mapped"

    log_time "Deinterleaving unmapped reads..."
    reformat.sh \
        in="$out_unmapped_raw".fq.gz \
        out1="$R1_unmapped" out2="$R2_unmapped"
    echo
else
    # Just move the files, if wanting to keep them interleaved
    mv -v "$out_mapped_raw".fq.gz "$outdir"/mapped
    mv -v "$out_unmapped_raw".fq.gz "$outdir"/unmapped
fi

# Move log files to main dir, remove temp files
log_time "Removing temporary files..."
mv "$outdir_sample"/mapped_raw/"$sample_id"*log "$LOG_DIR"
rm -rv "$outdir_sample"/mapped_raw "$outdir_sample"/unmapped_raw

# Quantify mapping success
n_mapped=$(zcat "$R1_mapped" | awk '{ s++ } END{ print s/4 }')
n_unmapped=$(zcat "$R1_unmapped" | awk '{ s++ } END{ print s/4 }')
pct=$(python3 -c "print(round($n_mapped / ($n_unmapped + $n_mapped) * 100, 2))")
log_time "Number of reads mapped/unmapped, and % mapped:\t$sample_id\t$n_mapped\t$n_unmapped\t$pct"

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
[[ "$deinterleave" = true ]] && ls -lh "$R1_mapped" "$R2_mapped" "$R1_unmapped" "$R2_unmapped"
final_reporting "$LOG_DIR"
