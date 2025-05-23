#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=star_index
#SBATCH --output=slurm-star_index-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Index a genome or transcriptome with STAR"
SCRIPT_VERSION="2023-08-13"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=STAR
TOOL_NAME=STAR
TOOL_DOCS="https://github.com/alexdobin/STAR, https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env_type=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/project/PAS0471/jelmer/conda/star
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
strict_bash=true

# Defaults - tool parameters
index_size="auto"
mem_bytes=4000000000

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
    echo "      sbatch $0 -i data/ref/genome.fa --annot data/ref/annotation.gtf -o results/star_index"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--infile         <file>  Input file: a nucleotide FASTA file with a genome or transcriptome assembly"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --annot             <file>  Reference annotation (GFF/GFF3/GTF) file (GTF preferred)  [default: no annotation file, but this is not recommended]"
    echo "  --index_size        <int>   Index size                              [default: $index_size => automatically determined from genome size]"
    echo "  --read_len          <int>   Read length (only applies with --annot) [default: unset]"
    echo "                              This will determine the overhang length, which is by default 100-1 = 99 bp."
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env_type               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env_type]"
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
    echo "NOTES:"
    echo "  The script will check how much memory has been allocated to the SLURM job (default: 64GB),"
    echo "  and pass that to STAR via the 'limitGenomeGenerateRAM argument'."
    echo "  When allocating more memory to the SLURM job,"
    echo "  wich can be necessary for large genomes, this will therefore be passed to STAR as well."
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
infile=
annot= && annot_opt=
read_len= && overhang_opt=
outdir=
opts=
version_only=false
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --infile )     shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --annot )           shift && annot=$1 ;;
        --index_size )      shift && index_size=$1 ;;
        --read_len )        shift && read_len=$1 ;;
        --opts )            shift && opts=$1 ;;
        --env_type )             shift && env_type=$1 ;;
        --no_strict )       strict_bash=false ;;
        --dl_container )    dl_container=true ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 && dl_container=true ;;
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
[[ "$strict_bash" == true ]] && set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input file specified, do so with -i/--infile" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input file $infile does not exist"
[[ -n "$annot" && ! -f "$annot" ]] && die "Annotation file $annot does not exist" "$all_opts"

# Define outputs and final ops based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
[[ "$IS_SLURM" == true ]] && mem_bytes=$((SLURM_MEM_PER_NODE * 1000000))
[[ -n "$annot" ]] && annot_opt="--sjdbGTFfile $annot"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Input assembly FASTA:                     $infile"
echo "Output dir:                               $outdir"
[[ -n "$annot" ]] && echo "Input annotation file:                    $annot"
[[ -n "$read_len" ]] && echo "Read length (for overhang size):                              $read_len"
[[ "$index_size" != "auto" ]] && echo "Index size:                               $index_size"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:        $opts"
log_time "Listing the input file(s):"
ls -lh "$infile"
[[ -n "$annot" ]] && ls -lh "$annot"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# STAR doesn't accept zipped FASTA files -- unzip if needed
if [[ $infile = *gz ]]; then
    infile_unzip=${infile/.gz/}
    if [[ ! -f $infile_unzip ]]; then
        log_time "Unzipping the currently gzipped FASTA file..."
        gunzip -c "$infile" > "$infile_unzip"
    else
        log_time "Using unzipped version of the FASTA file"
        ls -lh "$infile_unzip"
    fi
    infile="$infile_unzip"
fi

# Determine index size
if [[ "$index_size" == "auto" ]]; then
    log_time "Automatically determining the index size..."
    genome_size=$(grep -v "^>" "$infile" | wc -c)
    index_size=$(python -c "import math; print(math.floor(math.log($genome_size, 2)/2 -1))")
    log_time "Genome size (autom. determined):  $genome_size"
    log_time "Index size (autom. determined):   $index_size"
fi

# If a GTF/GFF file is provided, build the appropriate argument for STAR
if [[ -n "$read_len" ]]; then
    # Overhang length should be read length minus 1 - only if annot is included
    overhang=$(( read_len - 1 ))
    overhang_opt="--sjdbOverhang $overhang"
    log_time "Based on read length $read_len, setting overhang to: $overhang"
fi

log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM "$mem_bytes" \
    --genomeDir "$outdir" \
    --genomeFastaFiles "$infile" \
    --genomeSAindexNbases "$index_size" \
    --runThreadN "$threads" \
    $annot_opt \
    $overhang_opt \
    $opts

#? 2023-08-13: Removed this option: '--sjdbGTFtagExonParentTranscript Parent'

log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
final_reporting "$LOG_DIR"
