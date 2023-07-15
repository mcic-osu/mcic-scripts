#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=dl_genomes
#SBATCH --output=slurm-dl_genomes-%j.out

# Download genomes (and associated proteomes, annotations, etc) with the NCBI datasets tool

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants
readonly SCRIPT_NAME=dl-genomes.sh
readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="Jelmer Poelstra"
readonly SCRIPT_URL=https://github.com/mcic-osu/mcic-scripts
readonly MODULE=miniconda3/4.12.0-py39
readonly CONDA_ENV=/fs/ess/PAS0471/jelmer/conda/ncbi-datasets
readonly TOOL_BINARY="datasets download genome"
readonly TOOL_NAME="NCBI datasets"
readonly TOOL_DOCS=https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/

# Option defaults
include=genome
ref_only=false && ref_arg=
assembly_version=latest
is_slurm=true

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo
    echo "        $0 (v. $SCRIPT_VERSION): Run $TOOL_NAME"
    echo "        =============================================="
    echo "DESCRIPTION:"
    echo "  Download genomes (and associated proteomes, annotations, etc) with the NCBI datasets tool"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Download all C. elegans genomes:"
    echo "      sbatch $0 --taxon 'Caenorhabditis elegans' -o results/refgenome"
    echo "  - Download all reference nematode genomes:"
    echo "      sbatch $0 --taxon 'nematoda' --reference -o results/refgenomes"
    echo "  - Download a list of accessions:"
    echo "      sbatch $0 --accessions metadata/accession.txt -o results/refgenomes"
    echo "  - Also download annotation and proteome:"
    echo "      sbatch $0 --taxon 'human' --include 'genome,protein,gff' -o results"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  To specify genomes to download, use one of the following two options:"
    echo "  1) -a/--accessions   <file>  Text file with list of NCBI accession numbers, one per line"
    echo "  2) -t/--taxon        <str>   Taxon string, e.g. 'bos taurus' or 'nematoda'"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --include           <str>   Comma-separated string with one or more of the following options: [default: 'genome']"
    echo "                              'genome', 'rna', 'protein', 'cds', 'gff3', 'gtf', 'gbff', 'seq-report'"
    echo "  --assembly_version  <str>   'latest' or 'all'                       [default: 'lastest']"
    echo "  --ref_only                  When specifying a taxon, only download 'reference genomes' [default: download all matching genomes]"
    echo "  --more_args         <str>   Quoted string with more argument(s) for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help               Print this help message and exit"
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
    export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
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
tool_version() {
    set +e
    load_tool_conda
    $TOOL_BINARY --version
    set -e
}

# Print the tool's help
tool_help() {
    load_tool_conda
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
accessions= && accession_arg=
taxon= && taxon_arg=
outdir=
more_args=

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )         shift && outdir=$1 ;;
        -t | --taxon )          shift && taxon=$1 ;;
        -a | --accessions )     shift && accessions=$1 ;;
        --include )             shift && include=$1 ;;
        --assembly_version )    shift && assembly_version=$1 ;;
        --ref_only )            ref_only=true ;;
        --more_args )           shift && more_args=$1 ;;
        -v )                    script_version; exit 0 ;;
        -h | --help )           script_help; exit 0 ;;
        --version )             tool_version; exit 0 ;;
        * )                     die "Invalid option $1" "$all_args" ;;
    esac
    shift
done

# Check input
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_args"

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Check if this is a SLURM job
if [[ -z "$SLURM_JOB_ID" ]]; then is_slurm=false; else is_slurm=true; fi

# Strict bash settings
set -euo pipefail

# Load software
load_tool_conda

# ==============================================================================
#              DEFINE OUTPUTS AND DERIVED INPUTS, BUILD ARGS
# ==============================================================================
# Make paths absolute
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"

# Define outputs based on script parameters
readonly version_file="$outdir"/logs/version.txt
readonly log_dir="$outdir"/logs
[[ "$ref_only" == true ]] && ref_arg="--reference"

if [[ -n "$accessions" ]]; then
    accessions=$(realpath "$accessions")
    [[ ! -f "$accessions" ]] && die "Input file $accessions does not exist"
    subcommand=accession
    accession_arg="--inputfile $accessions"
elif [[ -n "$taxon" ]]; then
    subcommand=taxon
    taxon_arg="$taxon"
else
    die "Specify either an accession file or taxon string" "$all_args"
fi

# ==============================================================================
#                               REPORT
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Output dir:                               $outdir"
echo "Subcommand:                               $subcommand"
[[ -n "$taxon" ]] && echo "Taxon:                                    $taxon"
[[ -n "$accessions" ]] && echo "Input file with accession list:           $accessions"
[[ -n "$accessions" ]] && echo "Number of input accessions:               $(wc -l < "$accessions")"
[[ -n $more_args ]] && echo "Other arguments for $TOOL_NAME:   $more_args"
echo
[[ -n "$accessions" ]] && log_time "Listing the input file(s):"
[[ -n "$accessions" ]] && ls -lh "$accessions"
[[ "$is_slurm" = true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directories
log_time "Creating the output directories..."
mkdir -pv "$log_dir"

# Move into the output dir
cd "$outdir" || exit 1

# Report info on genomes #TODO
#info_fields="accession,assminfo-name,organism-name,assminfo-refseq-category,assminfo-level,assmstats-number-of-contigs,assmstats-contig-n50"
#log_time datasets summary genome $subcommand "$taxon_arg" $accession_arg --as-json-lines |
#    dataformat tsv genome --fields "$info_fields" |
#    column -ts $'\t'

# Run the tool
log_time "Running $TOOL_NAME..."
runstats $TOOL_BINARY $subcommand \
        ${accession_arg}"${taxon_arg}" \
        --include "$include" \
        --assembly-version "$assembly_version" \
        $ref_arg \
        --filename genomes.zip \
        --no-progressbar \
        --api-key "$NCBI_API_KEY" \
        $more_args

# Process the output - unzip
log_time "Unzipping the downloaded archive..."
unzip -o genomes.zip

# Process the output - move the genomic data files
if [[ $include =~ "genome" ]]; then
    log_time "Moving genome nucleotide FASTA files..."
    find ncbi_dataset -type f -name "*genomic.fna" -exec mv -v {} . \;
    log_time "Number of output genomes ('*fna' files):      $(find . -name "*genomic.fna" | wc -l)"
fi
if [[ $include =~ "protein" ]]; then
    log_time "Moving proteome FASTA files..."
    for proteome_old in $(find ncbi_dataset -type f -name "protein.faa"); do
        acc_nr=$(dirname "$proteome_old" | xargs basename)
        mv -v "$proteome_old" "$acc_nr"_protein.faa
    done
    log_time "Number of output proteomes:                   $(find . -name "*_protein.faa" | wc -l)"
fi
if [[ $include =~ "gff3" ]]; then
    log_time "Moving GFF files..."
    for proteome_old in $(find ncbi_dataset -type f -name "genomic.gff"); do
        acc_nr=$(dirname "$proteome_old" | xargs basename)
        mv -v "$proteome_old" "$acc_nr".gff
    done
    log_time "Number of output GFF files:                   $(find . -name "*.gff" | wc -l)"
fi
if [[ $include =~ "gtf" ]]; then
    log_time "Moving GTF files..."
    for proteome_old in $(find ncbi_dataset -type f -name "genomic.gtf"); do
        acc_nr=$(dirname "$proteome_old" | xargs basename)
        mv -v "$proteome_old" "$acc_nr".gtf
    done
    log_time "Number of output GTF files:                   $(find . -name "*.gtf" | wc -l)"
fi

# Report
[[ -n "$accessions" ]] && log_time "Number of input accessions:                   $(wc -l < "$accessions")"

# Clean up
log_time "Removing original ZIP file and NCBI's README file..."
rm -v README.md genomes.zip

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

#? Alternative:
# Use AstroBioMike's bit - https://github.com/AstrobioMike/bit
#micromamba activate /fs/ess/PAS0471/jelmer/conda/bit
#bit-dl-ncbi-assemblies -w "$acc_file" -f fasta -j 1
