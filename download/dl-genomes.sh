#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=dl_genomes
#SBATCH --output=slurm-dl_genomes-%j.out

#TODO - Implement '--prefer_refseq' option where RefSeq genomes will be downloaded
#TODO   as available, and GenBank genomes only for those for which there is no RefSeq genome

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Download genomes (and associated proteomes, annotations, etc) with the NCBI datasets tool"
SCRIPT_VERSION="2023-08-24"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh
TOOL_BINARY="datasets download genome"
TOOL_NAME="NCBI datasets"
TOOL_DOCS=https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets
VERSION_COMMAND="$TOOL_BINARY --version"

# Defaults - generics
env=conda                           # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/ncbi-datasets
container_path=
container_url=
dl_container=false
container_dir="$HOME/containers"
version_only=false

# Defaults - tool parameters
include="all"
ref_only=false
meta_only=false
assembly_version=latest
assembly_source=GenBank
move_output=true
meta_fields="accession,assminfo-name,organism-name,assminfo-refseq-category,assminfo-level,assmstats-number-of-contigs,assmstats-contig-n50,assminfo-sequencing-tech"

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
script_help() {
    echo "                          $0"
    echo "      (v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL)"
    echo "        =============================================================="
    echo "DESCRIPTION:"
    echo "  $DESCRIPTION"
    echo
    echo "USAGE / EXAMPLE COMMANDS:"
    echo "  - Download all C. elegans genomes:"
    echo "      sbatch $0 --taxon 'Caenorhabditis elegans' -o results/refgenome"
    echo "  - Download all reference nematode genomes:"
    echo "      sbatch $0 --taxon 'nematoda' --ref_only -o results/refgenomes"
    echo "  - Download a list of accessions:"
    echo "      sbatch $0 --accession_file metadata/accession.txt -o results/refgenomes"
    echo "  - Also download annotation and proteome:"
    echo "      sbatch $0 --taxon 'human' --include 'genome,protein,gff3' -o results"
    echo "  - Download a single accession:"
    echo "      sbatch $0 --accession GCA_000001405.29 -o data/ref --assembly_source all"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "To specify genomes to download, use _one_ of the following three options:"
    echo "  -t/--taxon          <str>   Taxon string (e.g. 'bos taurus' or 'nematoda') or NCBI Taxonomy ID (e.g. '10116')"
    echo "  -a/--accession_file <file>  Text file with list of NCBI accession numbers, one per line"
    echo "  -A/--accession      <str>   A single accession number to download (e.g. 'GCA_003693625.1')"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --as_is                     Don't move all files into the output dir - keep subdir structure with one folder per genome"         
    echo "  --include           <str>   Comma-separated string with one or more of the following options: [default: $include]"
    echo "                              Options: 'all' (all of the following), 'genome', 'rna', 'protein', 'cds', 'gff3', 'gtf', 'gbff', 'seq-report'"
    echo "  --assembly_source   <str>   'GenBank', 'RefSeq', or 'all'           [default: $assembly_source]"
    echo "                                - 'GenBank' will likely return more records"
    echo "                                - 'RefSeq' will likely return more annotations and proteomes" 
    echo "                                - Does not apply when you're directly requestion GCA_ (GenBank) or GCF_ (RefSeq) accession nrs"
    echo "  --assembly_version  <str>   'latest' or 'all'                       [default: $assembly_version]"
    echo "  --ref_only                  When specifying a taxon, only download 'reference genomes' [default: download all matching genomes]"
    echo "  --meta_fields       <str>   Comma-separated list of metadata fields for the 'selected' metadata file"
    echo "                                (Note: a metadata file with all possible metadata will also be produced.)"
    echo "                                [default: $meta_fields]"
    echo "  --meta_only                 Don't download genomes, only metadata   [default: $meta_only]"
    echo "  --opts              <str>   Quoted string with additional options for $TOOL_NAME"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --env               <str>   Use a Singularity container ('container') or a Conda env ('conda') [default: $env]"
    echo "  --container_url     <str>   URL to download the container from      [default: $container_url]"
    echo "                                A container will only be downloaded if an URL is provided with this option, or '--dl_container' is used"
    echo "  --container_dir     <str>   Dir to download the container to        [default: $container_dir]"
    echo "  --dl_container              Force a redownload of the container     [default: $dl_container]"
    echo "  --conda_env         <dir>   Full path to a Conda environment to use [default: $conda_path]"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
    echo "  --version                   Print the version of $TOOL_NAME and exit"
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
outdir=
accession_file=
accession=
taxon=
opts=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -o | --outdir )         shift && outdir=$1 ;;
        -t | --taxon )          shift && taxon=$1 ;;
        -a | --accession_file ) shift && accession_file=$1 ;;
        -A | --accession )      shift && accession=$1 ;;
        --include )             shift && include=$1 ;;
        --assembly_version )    shift && assembly_version=$1 ;;
        --assembly_source )     shift && assembly_source=$1 ;;
        --ref_only )            ref_only=true ;;
        --meta_fields )         shift && meta_fields=$1 ;;
        --as_is )               move_output=false ;;
        --meta_only )           meta_only=true ;;
        --opts )                shift && opts=$1 ;;
        --env )                 shift && env=$1 ;;
        --dl_container )        dl_container=true ;;
        --container_dir )       shift && container_dir=$1 ;;
        --container_url )       shift && container_url=$1 && dl_container=true ;;
        -h | --help )           script_help; exit 0 ;;
        -v )                    script_version; exit 0 ;;
        --version )             version_only=true ;;
        * )                     die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Load software
load_env "$conda_path" "$container_path" "$dl_container"
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
[[ "$version_only" == true ]] && tool_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"

# Define outputs based on script parameters
[[ ! "$outdir" =~ ^/ ]] && outdir="$PWD"/"$outdir"       # Make absolute because we will move into the outdir
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
meta_dir="$outdir"/metadata && mkdir -p "$meta_dir"
meta_all="$meta_dir"/meta_all.tsv
meta_sel="$meta_dir"/meta_sel.tsv

# Which seq files to include
[[ "$include" == "all" ]] && include="cds,gbff,genome,gff3,gtf,protein,rna,seq-report"

# Build command to specify which genomes to download
if [[ -n "$accession_file" ]]; then
    accession_file=$(realpath "$accession_file")
    [[ ! -f "$accession_file" ]] && die "Input file $accession_file does not exist"
    data_arg=(accession --inputfile "$accession_file")
elif [[ -n "$accession" ]]; then
    data_arg=(accession "$accession")
    # If accession nr is GCA_ or GCF_, set --assembly-source to all
    if echo "$accession" | grep -q "GC[AF]_"; then
        log_time "The accession nr. contains GCA_/GCF_, so setting --assembly-source to 'all'"
        assembly_source=all
    fi
elif [[ -n "$taxon" ]]; then
    data_arg=(taxon "$taxon")
else
    die "Use either --accession, --accession_file, or --taxon to specify which genomes to download" "$all_opts"
fi

# Complete data arg
[[ "$ref_only" == true ]] && data_arg+=("--reference")
data_arg+=(--assembly-version "$assembly_version")
data_arg+=(--assembly-source "$assembly_source")

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Output dir:                               $outdir"
[[ -n "$taxon" ]] && echo "Taxon:                                    $taxon"
[[ -n "$accession_file" ]] && echo "Input file with accession list:           $accession_file"
[[ -n "$accession_file" ]] && echo "Number of input accessions:               $(wc -l < "$accession_file")"
echo "What to download for each genome:         $include"
echo "Assembly source (GenBank vs RefSeq):      $assembly_source"
echo "Download 'reference' genomes only:        $ref_only"
echo "Which assembly version(s) to download:    $assembly_version"
echo "Download metadata only?                   $meta_only"
[[ -n $opts ]] && echo "Additional options for $TOOL_NAME:        $opts"
[[ -n "$accession_file" ]] && log_time "Listing the input file(s):"
[[ -n "$accession_file" ]] && ls -lh "$accession_file"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Move into the output dir
cd "$outdir" || exit 1

# Get metadata on genomes
log_time "Storing *full* metadata on genomes in $meta_all..."
runstats datasets summary genome "${data_arg[@]}" --as-json-lines |
    dataformat tsv genome > "$meta_all"

log_time "Storing *selected* metadata on genomes in $meta_sel..."
runstats datasets summary genome "${data_arg[@]}" --as-json-lines |
    dataformat tsv genome --fields "$meta_fields" > "$meta_sel"

# Report metadata
n_genomes=$(tail -n +2 "$meta_sel" | wc -l)
log_time "Found $n_genomes genomes"
log_time "Showing genome metadata file $meta_sel"
column -ts $'\t' "$meta_sel"

# Download the genomes
if [[ "$meta_only" == false ]]; then
    log_time "Downloading the genomes..."
    runstats $CONTAINER_PREFIX $TOOL_BINARY \
        "${data_arg[@]}" \
        --include "$include" \
        --filename genomes.zip \
        --no-progressbar \
        --api-key "$NCBI_API_KEY" \
        $opts

    # Process the output - unzip the archive
    log_time "Unzipping the downloaded archive..."
    unzip -q -o genomes.zip

    # Process the output - move the genomic data files
    if [[ "$move_output" == true ]]; then
        log_time "Moving and renaming the downloaded files..."

        while IFS= read -r -d '' file; do
            acc_nr=$(dirname "$file" | xargs basename)
            if [[ $(basename "$file") == "cds_from_genomic.fna" ]]; then
                # Make sure cds_from_genomic.fna is distinct from the genome fna
                outfile="$acc_nr"_cds.fasta
            elif [[ $(basename "$file") == "rna.fna" ]]; then
                outfile="$acc_nr"_rna.fasta
            else
                outfile="$acc_nr"."${file##*.}"
            fi
            mv -v "$file" "$outfile"
        done < <(find ncbi_dataset -type f -wholename "*data/GC*" -print0)

        echo
        echo "# Nr output nucleotide FASTAs:          $(find . -name "*.fna" | wc -l)"
        echo "# Nr output proteomes:                  $(find . -name "*.faa" | wc -l)"
        echo "# Nr output GFF files:                  $(find . -name "*.gff" | wc -l)"
        echo "# Nr output GTF files:                  $(find . -name "*.gtf" | wc -l)"
    fi

    # Report nr of input accessions again, for comparison
    [[ -n "$accession_file" ]] && log_time "Number of input accessions:                   $(wc -l < "$accession_file")"

    # Clean up
    log_time "Removing original ZIP file and NCBI's README file..."
    rm -r README.md genomes.zip
fi

# Show output metadata files
log_time "Output metadata files:"
ls -lh "$meta_all" "$meta_sel"

final_reporting "$LOG_DIR"
