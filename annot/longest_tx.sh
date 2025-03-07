#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=longest_tx
#SBATCH --output=slurm-longest_tx-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Create a longest-transcript/isoform proteome FASTA
Works with matching NCBI proteome FASTA and GTF files
Also works with proteome FASTA files only, as long as transcripts/isoforms are
IDed by .t1/.p1-style suffices"
SCRIPT_VERSION="2023-12-18"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_NAME=seqkit
VERSION_COMMAND="seqkit | sed -n '3p'"

# Defaults - generics
conda_path=/fs/ess/PAS0471/jelmer/conda/seqkit
strict_bash=true
version_only=false                  # When true, just print tool & script version info and exit

# Defaults - tool parameters
keep_intermed=false                 # Remove intermediate files
t1_style=true                       # Isoforms have .t1/.p1 suffices in proteome file
                                    # Currently the only supported format

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
    echo "  - Basic usage:"
    echo "      sbatch $0 -i my.faa --gtf my.gtf -o results/longest_iso"
    echo "      sbatch $0 -i my.faa -o results/longest_iso"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--faa            <file>  Input proteome FASTA file with multiple isoforms per gene"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --gtf               <file>  Input annotation file in GTF format"
    echo "  --keep_intermed             Keep intermediate files [default: remove]"
    echo
    echo "UTILITY OPTIONS:"
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
infile=
gtf=
outdir=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --faa )        shift && infile=$1 ;;
        --gtf )             shift && gtf=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        --keep_intermed )   keep_intermed=true ;;
        --no_strict )       strict_bash=false ;;
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
#                               RUN
# ==============================================================================
# Clean the FASTA headers
#if [[ "$clean_header" == true ]]; then
#    log_time "Keeping only the first 'word' in the FASTA header..."
#    awk '{print $1}' "$infile" > "$infile_prepped"
#else
#    infile_prepped="$infile"
#fi
#else
#    log_time "Replacing tabs with spaces in the FASTA header..."
#    tr "\t" " " < "$infile" > "$infile_prepped"
#fi

log_time "Getting the isoform lengths..."
awk '{print $1}' "$infile" | seqkit fx2tab --length --name | sort -k1,1 > "$iso_lens"

log_time "Creating a gene2isoform lookup file..."
if [[ -n "$gtf" ]]; then
    # Use a GTF file to ID genes and proteins #TODO also allow for transcript_id
    grep -v "^#" "$gtf" | awk '$3 == "CDS"' |
        sed -E 's/.*gene_id "([^;]+)";.*protein_id "([^;]+)";.*/\1\t\2/' |
        sort -u | sort -k2,2 > "$gene2iso"
    n_genes_gtf=$(grep -v "^#" "$gtf" | awk '$3 == "gene"' | wc -l)
    log_time "Number of genes in the GTF file (can include noncoding): $n_genes_gtf"
elif [[ "$t1_style" == true ]]; then
    # Extract from the FASTA file directly
    awk -v OFS="\t" '{print $1,$1}' "$iso_lens" | sed -E 's/\.[pt][0-9]+//' |
        sort -k2,2 > "$gene2iso"
fi

n_genes=$(cut -f1 "$gene2iso" | sort -u | wc -l)
log_time "Number of genes in the original proteome file: $n_genes"

log_time "Getting the IDs of the longest isoforms..."
join -t$'\t' -1 2 -2 1 "$gene2iso" "$iso_lens" |
    sort -k2,2 -u | cut -f1 > "$longest_iso_ids"

log_time "Extracting the longest isoforms..."
seqkit grep -f "$longest_iso_ids" "$infile" > "$outfile"

# Report
n_in=$(grep -c "^>" "$infile")
n_ids=$(wc -l < "$longest_iso_ids")
n_out=$(grep -c "^>" "$outfile")

log_time "Number of proteins in the input file:         $n_in"
log_time "Number of proteins in the output file:        $n_out"

# Check that all proteins were found in the input file
if [[ ! "$n_ids" -eq "$n_out" ]]; then
    log_time "WARNING: Nr of longest-isoform IDs ($n_ids) is not the same as the
    nr of proteins in the output file ($n_out)"
fi

# Remove intermediate files
if [[ "$keep_intermed" == false ]]; then
    rm "$gene2iso" "$iso_lens" "$longest_iso_ids"
fi

log_time "Listing the output file:"
ls -lh "$outfile"
final_reporting "$LOG_DIR"
