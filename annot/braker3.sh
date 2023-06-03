#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=braker3
#SBATCH --output=slurm-braker3-%j.out

# Run Braker3 to annotate a genome assembly
# As 'evidence' for annotation, provide RNAseq data of the same species and/or
# protein data of the same or related species 

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=braker3.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly TOOL_BINARY_PART=braker.pl
readonly TOOL_NAME=Braker3
readonly TOOL_DOCS=https://github.com/Gaius-Augustus/BRAKER
readonly TOOL_PAPER=https://academic.oup.com/nargab/article/3/1/lqaa108/6066535

# Option defaults
container_path=/fs/ess/PAS0471/containers/braker3.sif
rebuild_container=false

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "======================================================================"
    echo "DESCRIPTION:"
    echo "  Run Braker3 to annotate a genome, using RNAseq and/or protein data"
    echo 
    echo "USAGE:"
    echo "  sbatch $0 -i <input FASTA> -o <outdir> --species <species_name> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--assembly       <file>  Input assembly nucleotide FASTA"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  --species           <str>   Species name (without space, e.g. 'homo_sapiens')"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --rnaseq_bamdir     <dir>   Directory/ies with RNAseq FASTQ or BAM files"
    echo "  --prot_seq          <file>  FASTA with reference proteins. For info on how to create this file:"
    echo "                              https://github.com/gatech-genemark/ProtHint#protein-database-preparation"
    echo " --rebuild_container          Rebuild the container from docker://teambraker/braker3:latest, e.g. to update to latest version"
    echo " --container_path     <file>  Use the specified container SIF file    [default: use /fs/ess/PAS0471/containers/braker3.sif]"
    echo "  --more_args         <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h                          Print this help message and exit"
    echo "  --help                      Print the help for $TOOL_NAME and exit"
    echo "  -v                          rint the version of this script and exit"
    echo "  -v/--version                Print the version of $TOOL_NAME and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i results/genome.fa -o results/braker --prot_seq data/ref/proteins.faa --species homo_sapiens"
    echo "  sbatch $0 -i results/genome.fa -o results/braker --prot_seq data/ref/proteins.faa --species homo_sapiens --more_args '--fungus'"
    echo
    echo "HARDCODED PARAMETERS:"
    echo "  - '--verbosity=3'           Set Braker verbosity level to 3"
    echo
    echo "TOOL DOCUMENTATION:"
    echo "  - Docs: $TOOL_DOCS"
    echo "  - Paper: $TOOL_PAPER"
    echo
}

load_tool() {
    # GeneMark license key
    # See https://github.com/Gaius-Augustus/BRAKER#genemark-ex
    GENEMARK_BASEDIR=/fs/project/PAS0471/jelmer/software/genemark-ex
    cp "$GENEMARK_BASEDIR"/gm_key_64 ~/.gm_key
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
tool_version() {
    $TOOL_BINARY --version
}

# Print the tool's help
tool_help() {
    $TOOL_BINARY --help
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
assembly=""
species=""
outdir=""
prot_seq="" && prot_seq_arg=""
rnaseq_bamdir="" && rna_seq_arg=""
more_args=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --assembly )       shift && assembly=$1 ;;
        -o | --outdir )         shift && outdir=$1 ;;
        --species )             shift && species=$1 ;;
        --prot_seq )            shift && prot_seq=$1 ;;
        --rnaseq_bamdir )       shift && rnaseq_bamdir=$1 ;;
        --rebuild_container )   readonly rebuild_container=true ;;
        --container_path )      readonly container_path=$1 ;;
        --more_args )           shift && readonly more_args=$1 ;;
        -v )                    script_version; exit 0 ;;
        -h )                    script_help; exit 0 ;;
        --version )             tool_version; exit 0 ;;
        --help )                tool_help; exit 0;;
        * )                     die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$assembly" ]] && die "No input assembly specified, do so with -i/--assembly" "$all_args"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"
[[ -z "$species" ]] && die "Please specify a species name with --species" "$all_args"
[[ ! -f "$assembly" ]] && die "Input assembly $assembly does not exist"
[[ -n "$prot_seq" && ! -f "$prot_seq" ]] && die "Protein seq file $prot_seq does not exist" 

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Load software and set nr of threads
set_threads

# Make input paths absolute
[[ ! $assembly =~ ^/ ]] && assembly=$(realpath "$assembly")
[[ -n "$prot_seq" ]] && [[ ! $prot_seq =~ ^/ ]] && prot_seq=$(realpath "$prot_seq")
[[ -n "$rnaseq_bamdir" ]] && [[ ! $rnaseq_bamdir =~ ^/ ]] && rnaseq_bamdir=$(realpath "$rnaseq_bamdir")

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Protein sequence file
[[ -n "$prot_seq" ]] && prot_seq_arg="--prot_seq=$prot_seq"

# RNAseq BAM files
if [[ -n "$rnaseq_bamdir" ]]; then
    bam_list=$(find "$rnaseq_bamdir" -name "*bam" | tr "\n" "," | sed 's/,$/\n/')
    rna_seq_arg="--bam=$bam_list"
fi

# Augustus config dir
augustus_config_dir="$PWD"/"$outdir"/augustus_config
augustus_config_arg="--AUGUSTUS_CONFIG_PATH=$augustus_config_dir"

# Define outputs based on script parameters
readonly TOOL_BINARY="singularity exec $container_path $TOOL_BINARY_PART"
readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Output dir:                               $outdir"
echo "Input file:                               $assembly"
echo "Species:                                  $species"
[[ -n "$prot_seq" ]] && echo "Reference protein FASTA:                  $prot_seq"
[[ -n "$rnaseq_bamdir" ]] && echo "RNAseq data dir:                          $rnaseq_bamdir"
[[ $more_args != "" ]] && echo "Other arguments for $TOOL_NAME:           $more_args"
echo "Number of threads/cores:                  $threads"
echo
echo "Listing the input file(s):"
ls -lh "$assembly"
[[ -n "$prot_seq" ]] && ls -lh "$prot_seq"
[[ -n "$rnaseq_bamdir" ]] && ls -lh "$rnaseq_bamdir"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
log_time "Creating the output directories..."
mkdir -pv "$log_dir"

# Move into the output dir
#log_time "Moving into the output dir $outdir..."
#cd "$outdir" || exit 1

# Rebuild the container SIF if needed
if [[ "$rebuild_container" == "true" ]]; then
    log_time "Rebuilding container SIF file from the internet..."
    singularity build "$container_path" docker://teambraker/braker3:latest
fi

# Copy the Augustus config dir from the container
# This is needed because Braker will otherwise try to write inside the container
log_time "Copying the Augustus config dir from the container..."
singularity exec -B "$PWD":"$PWD" "$container_path" \
    cp -r /usr/share/augustus/config/ "$augustus_config_dir"

# Run
log_time "Running $TOOL_NAME..."
runstats \
    $TOOL_BINARY \
    --workingdir="$outdir" \
    --genome="$assembly" \
    $prot_seq_arg \
    $rna_seq_arg \
    $augustus_config_arg \
    --species="$species" \
    --verbosity=3 \
    --threads="$threads" \
    $more_args

# Other useful options:
#  --gff3               Output a GFF3 instead of a GTF file
#  --fungus             GeneMark-EX option: run algorithm with branch point model
#                       (most useful for fungal genomes)
#  --useexisting        Use the present config and parameter files if they exist for 'species'; 
#                       will overwrite original parameters if BRAKER performs an AUGUSTUS training.

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
printf "\n======================================================================"
log_time "Versions used:"
tool_version | tee "$version_file"
script_version | tee -a "$version_file" 
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
[[ "$is_slurm" = true ]] && echo && resource_usage
log_time "Done with script $SCRIPT_NAME"
echo
