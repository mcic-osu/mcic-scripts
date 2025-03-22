#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --job-name=longest_tx
#SBATCH --output=slurm-longest_tx-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Create a longest-transcript/isoform proteome FASTA.
- Works with matching proteome FASTA and GTF files
- Also works with proteome FASTA files only, if transcripts/isoforms are
  IDed by .t1/.p1-style suffices
- If multiple isoforms have the same length, the first one is chosen
"
SCRIPT_VERSION="2025-03-22"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_NAME=seqkit
VERSION_COMMAND="seqkit | sed -n '3p'"

# Defaults - generics
env_type=conda
conda_path=/fs/ess/PAS0471/jelmer/conda/seqkit
container_dir="$HOME/containers"
container_url=
container_path=

# Defaults - tool parameters
keep_intermed=false                 # Remove intermediate files
t1_style=true                       # Isoforms have .t1/.p1 suffices in proteome file
                                    # Currently the only supported format

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "
                        $0
    v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL
            =================================================

DESCRIPTION:
$DESCRIPTION
    
USAGE / EXAMPLE COMMANDS:
  - Basic usage example:
      sbatch $0 -i my.faa --gtf my.gtf -o results/longest_iso
      sbatch $0 -i my.faa -o results/longest_iso
    
REQUIRED OPTIONS:
  -i/--faa            <file>  Input proteome FASTA file with multiple isoforms per gene
  -o/--outdir         <dir>   Output dir (will be created if needed)
    
OTHER KEY OPTIONS:
  --gtf               <file>  Input annotation file in GTF format
  --keep_intermed             Keep intermediate files                           [default: remove]
    
UTILITY OPTIONS:
  NOTE: The software used in this script is Seqkit (https://bioinf.shenwei.me/seqkit)
  --env_type          <str>   Use a Singularity container ('container')         [default: $env_type]
                              or a Conda environment ('conda') 
  --conda_path        <dir>   Full path to a Conda environment to use           [default: $conda_path]
  --container_dir     <str>   Dir to download a container to                    [default: $container_dir]
  --container_url     <str>   URL to download a container from                  [default (if any): $container_url]
  --container_path    <file>  Local singularity image file (.sif) to use        [default (if any): $container_path]
  -h/--help                   Print this help message
  -v/--version                Print script and $TOOL_NAME versions
"
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
        wget -q "$FUNCTION_SCRIPT_URL" -O "$function_script"
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
version_only=false                  # When true, just print tool & script version info and exit
infile=
gtf=
outdir=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --faa )        shift && infile=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --gtf )             shift && gtf=$1 ;;
        --keep_intermed )   keep_intermed=true ;;
        --env_type )        shift && env_type=$1 ;;
        --conda_path )      shift && conda_path=$1 ;;
        --container_dir )   shift && container_dir=$1 ;;
        --container_url )   shift && container_url=$1 ;;
        --container_path )  shift && container_path=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v | --version )    version_only=true ;;
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
load_env "$conda_path"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$infile" ]] && die "No input FASTA file specified, do so with -i/--faa" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -f "$infile" ]] && die "Input FASTA file $infile does not exist"
[[ -n "$gtf" && ! -f "$gtf" ]] && die "Input GTF file $gtf does not exist"

# Further processing
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
indir=$(dirname "$infile")
[[ "$indir" == "$outdir" ]] && die "Input dir can't be the same as the output dir ($indir)"

# Define outputs
file_id=$(basename "${infile%.*}")
gene2iso="$outdir"/"$file_id"_gene2iso.tsv
iso_lens="$outdir"/"$file_id"_iso_lens.tsv
longest_iso_ids="$outdir"/"$file_id"_longest_iso_ids.txt
outfile="$outdir"/"$file_id".faa

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Output dir:                               $outdir"
echo "Input FASTA file:                         $infile"
[[ -n "$gtf" ]] && echo "Input GTF file:                           $gtf"
log_time "Listing the input file(s):"
ls -lh "$infile" 
[[ -n "$gtf" ]] && ls -lh "$gtf" 
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               MAIN
# ==============================================================================
# Create a table with the length of each isoform
log_time "Getting the isoform lengths..."
awk '{print $1}' "$infile" |
    seqkit fx2tab --length --name |
    sort -k1,1 > "$iso_lens"

# Create a gene to isoform ID lookup table
log_time "Creating a gene2isoform lookup file..."
if [[ -n "$gtf" ]]; then
    # Use a GTF file to ID genes and proteins
    #TODO also allow for transcript_id
    grep -v "^#" "$gtf" |
        awk '$3 == "CDS"' |
        sed -E 's/.*gene_id "([^;]+)";.*protein_id "([^;]+)";.*/\1\t\2/' |
        sort -u |
        sort -k2,2 > "$gene2iso"
    
    n_genes_gtf=$(grep -v "^#" "$gtf" | awk '$3 == "gene"' | wc -l)
    log_time "Total number of genes in the GTF file quantified by counting
    'gene' entries in third column (can include noncoding): $n_genes_gtf"
elif [[ "$t1_style" == true ]]; then
    # Extract gene2iso lookup directly from the FASTA file
    awk -v OFS="\t" '{print $1,$1}' "$iso_lens" |
        sed -E 's/\.[pt][0-9]+//' |
        sort -k2,2 > "$gene2iso"
fi

# Get a list with the longest isoform for each gene
log_time "Getting the IDs of the longest isoforms..."
join -t$'\t' -1 2 -2 1 "$gene2iso" "$iso_lens" | 
    sort -k3,3nr |
    sort -k2,2 -u |
    cut -f1 > "$longest_iso_ids"

# Create the final output file, a protein FASTA containing only the longest
# isoform for each gene
log_time "Extracting the longest isoforms..."
seqkit grep -f "$longest_iso_ids" "$infile" > "$outfile"

# Report
n_iso=$(cut -f2 "$gene2iso" | sort -u | wc -l)
n_genes=$(cut -f1 "$gene2iso" | sort -u | wc -l)
n_ids=$(wc -l < "$longest_iso_ids")
n_in=$(grep -c "^>" "$infile")
n_out=$(grep -c "^>" "$outfile")
log_time "Done. Numbers of genes and isoforms:"
echo "Number of isoforms/genes:                      $n_iso / $n_genes"
echo "Number of entries in the in-/output file:      $n_in / $n_out"

# Check that all proteins were found in the input file,
# and that the number of isoforms is higher than the number of genes
if [[ ! "$n_ids" -eq "$n_out" ]]; then
    log_time "WARNING: Nr of longest-isoform IDs ($n_ids) is not the same as the
    nr of proteins in the output file ($n_out)"
fi
if [[ "$n_iso" -eq "$n_genes" ]]; then
    log_time "WARNING: The number of isoforms is the same as the number of genes"
fi

# Remove intermediate files
if [[ "$keep_intermed" == false ]]; then
    rm "$gene2iso" "$iso_lens" "$longest_iso_ids"
fi

# Final reporting
log_time "Listing the output file:"
ls -lh "$outfile"
final_reporting "$LOG_DIR"
