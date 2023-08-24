#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=localtree
#SBATCH --output=slurm-localtree-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Create a local phylogenetic tree with Parsnp"
SCRIPT_VERSION="2023-08-24"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions2.sh

# Defaults
color_column=pathovar               # Name of the metadata column to color tip labels by
tiplab_column=isolate               # Name of the metadata column with alternative tip labels

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
    echo "      sbatch $0 -i results/spades --ref_id UPB820 --bed my.bed -o results/localtree"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  -i/--indir          <dir>   Dir with FASTA files to align (should also contain the ref FASTA)"
    echo "  -r/--ref_id         <file>  Reference genome ID (filename without extension and dir)"
    echo "  --bed               <file>  BED file with genomic coordinates of focal region"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --root              <str>   ID of the genome to root the tree with  [default: no rerooting]"
    echo "  --meta              <file>  File with metadata to plot the tree"
    echo "  --color_column      <str>   Name of the metadata column to color tip labels by ('NULL' to ignore) [default: $color_column]"
    echo "  --tiplab_column     <str>   Name of the metadata column with alternative tip labels ('NULL' to ignore) [default: $tiplab_column]"
    echo "                                NOTE: --color_column and --tiplab_column are only used when a metadata file is provided"
    echo
    echo "UTILITY OPTIONS:"
    echo "  -h/--help                   Print this help message and exit"
    echo "  -v                          Print the version of this script and exit"
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
indir=
ref_id=
outdir=
meta= && meta_opt=
bedfile=
root= && root_opt=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -r | --ref_id )     shift && ref_id=$1 ;;
        --bed )             shift && bedfile=$1 ;;
        --meta )            shift && meta=$1 ;;
        --color_column )    shift && color_column=$1 ;;
        --tiplab_column )   shift && tiplab_column=$1 ;;
        --root )            shift && root=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v )                script_version; exit 0 ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Check options provided to the script
[[ -z "$indir" ]] && die "No input dir specified, do so with -i/--indir" "$all_opts"
[[ -z "$bedfile" ]] && die "No BED file specified, do so with --bed" "$all_opts"
[[ -z "$ref_id" ]] && die "No reference genome ID specified, do so with -r/--ref_id" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ ! -d "$indir" ]] && die "Input dir $indir does not exist"
[[ ! -f "$bedfile" ]] && die "BED file $bedfile does not exist"
[[ -n "$meta" && ! -f "$meta" ]] && die "Metadata file $meta does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"

# Inputs
mapfile -t genomes < <(find "$indir" -iname '*.fasta' -or -iname '*.fa' -or -iname '*.fna' | grep -v "$ref_id")
ref="$indir"/"$ref_id".fna
[[ ! -f "$ref" ]] && die "Input reference genome FASTA $ref does not exist"

# Outputs
locus_id=$(basename "$bedfile" .bed)
tree_org="$outdir"/parsnp.tree
tree=${tree_org/.tree/_fixnames.tree}

# Tree plot script options
[[ -n "$root" ]] && root_opt="--root $root"
[[ -n "$meta" ]] && meta_opt="--annot $meta --color_column $color_column --tiplab_column $tiplab_column"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Output dir:                               $outdir"
echo "Reference genome FASTA:                   $ref"
echo "Input dir with genomes:                   $indir"
echo "BED file with locus coordinates:          $bedfile"
[[ -n "$meta" ]] && echo "Metadata file:                            $meta"
[[ -n "$root" ]] && echo "Tree root:                                $root"
echo "Locus ID:                                 $locus_id"
echo "Number of input genomes:                  ${#genomes[@]}"
log_time "Listing the reference FASTA file:"
ls -lh "$ref"
log_time "Listing the input genomes:"
ls -lh "${genomes[@]}"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
module load miniconda3/23.3.1-py310

# Extract FASTA from BED file with all desired regions
log_time "Running bedtools getfasta..."
source activate /fs/ess/PAS0471/jelmer/conda/bedtools
bedtools getfasta -fi "$ref" -bed "$bedfile" -fo "$outdir"/"$locus_id".fa
log_time "Showing the bedtools output FASTA file:"
ls -lh "$outdir"/"$locus_id".fa

# Run Parsnp
log_time "Running Parsnp..."
conda deactivate && source activate /fs/ess/PAS0471/jelmer/conda/parsnp
runstats parsnp \
    --output-dir "$outdir" \
    --reference "$outdir"/"$locus_id".fa \
    --vcf \
    --curated \
    --threads "$threads" \
    --sequences "${genomes[@]}"
log_time "Showing the parsnp output tree file:"
ls -lh "$tree_org"

# Create a multi-FASTA alignment file (https://harvest.readthedocs.io/en/latest/content/harvest/quickstart.html)
log_time "Running Harvesttool to convert to a FASTA alignment file..."
runstats harvesttools -i "$outdir"/parsnp.ggr -M "$outdir"/parsnp.aln

# Fix the sample IDs in the tree, so they match the IDs in the metadata file
log_time "Fixing the sample IDs in the tree..."
sed -e 's/.fna//g' -e "s/$locus_id.fa.ref/$ref_id/g" -e "s/'//g" "$tree_org" > "$tree"
    
# Plot the tree
log_time "Plotting the tree..."
conda deactivate && source activate /fs/ess/PAS0471/jelmer/conda/r_tree
runstats Rscript mcic-scripts/trees/ggtree.R -i "$tree" $root_opt $meta_opt

# Report
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*
