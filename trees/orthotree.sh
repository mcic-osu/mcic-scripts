#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=orthotree
#SBATCH --output=slurm-orthotree-%j.out

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Create both a nucleotide and a protein phylogenetic tree for a gene
in a focal genome based on its Orthofinder orthogroup"
SCRIPT_VERSION="2023-10-18"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
conda_path=/fs/ess/PAS0471/jelmer/conda/iqtree          # Includes IQ-tree, MAFFT, csvtk

# Defaults
nboot=10000                          # Number of IQ-tree ultrafast bootstraps

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
    echo "      sbatch $0 -i results/genomes --ref_id UPB820 --bed my.bed --ortho_dir results/orthofinder -o results/orthotree"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--indir          <dir>   Dir with a FASTA (.fna) and GFF (.gff) file for each genome"
    echo "  -o/--outdir         <dir>   Output dir (will be created if needed)"
    echo "  --ortho_dir         <dir>   Dir with Orthofinder output for the same set of genomes as in the --indir"
    echo "  -r/--ref_id         <file>  Reference genome ID (filename without extension and dir)"
    echo "  --bed               <file>  BED file with genomic coordinates of focal region"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --nboot             <int>   Number of Ultrafast IQ-tree bootstraps  [default: $nboot]"
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
indir=
ref_id=
outdir=
bedfile=
boot_opt=
threads=

# Parse command-line args
all_opts="$*"
while [ "$1" != "" ]; do
    case "$1" in
        -i | --indir )      shift && indir=$1 ;;
        -r | --ref_id )     shift && ref_id=$1 ;;
        --bed )             shift && bedfile=$1 ;;
        --ortho_dir )       shift && ortho_dir=$1 ;;
        --nboot )           shift && nboot=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Strict Bash settings
set -euo pipefail

# Load conda env
load_env "$conda_path"

# Check options provided to the script
[[ -z "$indir" ]] && die "No input dir specified, do so with -i/--indir" "$all_opts"
[[ -z "$bedfile" ]] && die "No BED file specified, do so with --bed" "$all_opts"
[[ -z "$ref_id" ]] && die "No reference genome ID specified, do so with -r/--ref_id" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ -z "$ortho_dir" ]] && die "No Orthofinder output dir specified, do so with --ortho_dir" "$all_opts"
[[ ! -d "$indir" ]] && die "Input dir $indir does not exist"
[[ ! -d "$ortho_dir" ]] && die "Orthofinder output dir $ortho_dir does not exist"
[[ ! -f "$bedfile" ]] && die "BED file $bedfile does not exist"

# Define outputs based on script parameters
LOG_DIR="$outdir"/logs && mkdir -p "$LOG_DIR"
mem_gb=$((8*(SLURM_MEM_PER_NODE / 1000)/10))G   # 80% of available memory in GB

# Inputs
mapfile -t genomes < <(find "$indir" -iname '*.fasta' -or -iname '*.fa' -or -iname '*.fna' | grep -v "$ref_id")
ref_fna="$indir"/"$ref_id".fna
ref_gff_org="$indir"/"$ref_id".gff
orthogroups_file="$ortho_dir"/Orthogroups/Orthogroups.tsv
[[ ! -f "$ref_fna" ]] && die "Input reference genome FASTA $ref_fna does not exist"
[[ ! -f "$ref_gff_org" ]] && die "Input reference genome GFF $ref_gff_org does not exist"
[[ ! -f "$orthogroups_file" ]] && die "Orthogroups file $orthogroups_file does not exist"

# Outputs
locus_id=$(basename "$bedfile" .bed)
ref_gff="$outdir"/ref_gff_noseqs/$(basename "$ref_gff_org")
fna_unaln="$outdir"/fna_all/"$locus_id"_concat.fna
fna_aln="$outdir"/fna_all/"$locus_id"_aln.fna
asm_list="$outdir"/assembly_list.txt
msa="$outdir"/msa_prot.faa
root_file="$outdir"/orthofinder_root.txt

# Options
[[ "$nboot" -gt 0 ]] && boot_opt="--ufboot $nboot"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Output dir:                               $outdir"
echo "Reference genome FASTA:                   $ref_fna"
echo "Reference genome GFF:                     $ref_gff_org"
echo "Input dir with genomes:                   $indir"
echo "BED file with locus coordinates:          $bedfile"
echo "Orthofinder orthogroups file:             $orthogroups_file"
echo "Number of bootstraps:                     $nboot"
echo "Locus ID:                                 $locus_id"
echo "Number of input genomes:                  ${#genomes[@]}"
log_time "Listing the reference genomes files:"
ls -lh "$ref_fna" "$ref_gff_org"
log_time "Listing the first few input genomes:"
ls -lh "${genomes[@]:0:3}"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                           RUN - PREP
# ==============================================================================
# Create the output dirs
mkdir -p "$outdir"/fna_each "$outdir"/fna_all "$outdir"/trees "$outdir"/ref_gff_noseqs

# Make GFF without seqs for ref genome
log_time "Creating a GFF file for the reference without nucleotide sequences..."
sed '/^##FASTA/Q' "$ref_gff_org" > "$ref_gff"
ls -lh "$ref_gff"

# Get gene ID, orthogoup, and root # todo - check that there's exactly one match
log_time "Getting the focal gene ID and Orthogroup..."
gene_id=$(bedtools intersect -F 0.9 -a "$ref_gff" -b "$bedfile" | sed -E 's/.*ID=([^;]+);.*/\1/')

OG=$(grep "$gene_id" "$orthogroups_file" | cut -f1)
ortho_tree="$ortho_dir"/Gene_Trees/"$OG"_tree.txt

root=$(cut -d, -f1 "$ortho_tree" | cut -d: -f1 | sed -E 's/\(+//' | sed -E 's/_([0-9])_/.\1_/')
echo "$root" > "$root_file"
log_time "Gene ID: $gene_id   Orthogroup: $OG    Orthofinder root: $root"

# Create a table showing the number of gene copies in the focal orthogroup for each included genomes
OG_file1="$outdir"/OG.tsv               # File with all genes in the focal OG for each genome
OG_file2="$outdir"/OG_gene_nrs.tsv      # File with just the nr of genes in the focal OG in each genome

cat <(head -n1 "$orthogroups_file") <(grep "$OG" "$orthogroups_file") | \
    csvtk -t transpose | tail -n+2 | sed 's/, /\t/g' > "$OG_file1"
paste <(cut -f 1 "$OG_file1") <(awk '{print NF-1}' "$OG_file1") > "$OG_file2"

n1=$(awk '$2 == 1' "$OG_file2" | wc -l)
nmult=$(awk '$2 > 1' "$OG_file2" | wc -l)
nmiss=$(awk '$2 == 0' "$OG_file2" | wc -l)
log_time "Number of genomes with 1 / >1 / 0 copies: $n1 / $nmult / $nmiss"
log_time "See files $OG_file1 and $OG_file2"

# Assign the correct output tree file names, this depends on whether bootstraps will be done
if [[ "$nboot" -gt 0 ]]; then
    nuc_tree="$outdir"/trees/"$OG"_nuc.contree
    prot_tree="$outdir"/trees/"$OG"_prot.contree
else
    nuc_tree="$outdir"/trees/"$OG"_nuc.treefile
    prot_tree="$outdir"/trees/"$OG"_prot.treefile
fi

# ==============================================================================
#                           RUN - PROTEIN TREE
# ==============================================================================
echo -e "\n====================================================================="
msa_org="$ortho_dir"/MultipleSequenceAlignments/"$OG".fa
log_time "Listing the input MSA from Orthofinder, to build a protein tree:"
ls -lh "$msa_org"

log_time "Modifying the MSA IDs..."
sed -E '/>GC[AF]/s/_([0-9])_/.\1_/' "$msa_org" > "$msa"
ls -lh "$msa"

log_time "Running IQtree to build a protein tree..."
runstats iqtree \
    -s "$msa" --prefix "$outdir"/trees/"$OG"_prot $boot_opt \
    -nt "$threads" -ntmax "$threads" -mem "$mem_gb" -redo \
    > "$outdir"/logs/iqtree_"$OG"_prot.log
log_time "...Done. IQtree log:"
ls -lh "$outdir"/logs/iqtree_"$OG"_prot.log

# ==============================================================================
#                           RUN - NUCLEOTIDE TREE
# ==============================================================================
echo -e "\n====================================================================="
log_time "Looping over sequences in protein MSA file to create nucleotide FASTAs..."

>"$asm_list"

grep ">" "$msa" | while read -r seq_id; do
    # Get assembly and gene IDs
    asm_id=$(echo "$seq_id" | sed -E 's/>(.*)_([A-Z]+_[0-9]+)$/\1/' | sed -E 's/_([0-9])$/.\1/')
    gene_id=$(echo "$seq_id" | sed -E 's/>(.*)_([A-Z]+_[0-9]+)$/\2/')
    # Inputs
    gff_focal="$indir"/"$asm_id".gff
    asm_focal="$indir"/"$asm_id".fna
    # Outputs
    gene_fna="$outdir"/fna_each/"$asm_id"_"$gene_id".fna
    tmp_gff="$outdir"/fna_each/"$asm_id"_"$gene_id".tmp
    log_time "asm_id: $asm_id  /  gene_id: $gene_id"

    # Create nucleotide FASTA
    awk '$3 == "CDS"' "$gff_focal" | grep "$gene_id" > "$tmp_gff"
    if [[ -s "$tmp_gff" ]]; then
        # Extract FASTA by coords
        bedtools getfasta -fi "$asm_focal" -bed "$tmp_gff" > "$gene_fna"
        # Prepend assembly ID + gene ID to header
        sed -i "s/>/>${asm_id}_${gene_id} /" "$gene_fna"
        # Check nr of entries in the FASTA
        n_entry=$(grep -c "^>" "$gene_fna")
        [[ "$n_entry" -eq 0 ]] && echo "WARNING: FASTA $gene_fna is empty "
        if [[ "$n_entry" -gt 1 ]]; then
            echo "WARNING: FASTA $gene_fna contains >1 entries: $n_entry
            ls -lh $gene_fna" "$tmp_gff"
        else
            rm "$tmp_gff"
        fi
    else
        log_time "WARNING: Gene not found ($asm_id / $gene_id)"
    fi
    # Keep track of assemblies to check for duplicates
    echo "$asm_id" >> "$asm_list"
done

# Check for assemblies that are present twice
log_time "Printing assemblies with multiple copies of the focal gene (none if there's no output below):"
sort "$asm_list" | uniq -c | awk '$1 > 1' | sort -nr

# Concatenate FASTA files
log_time "Creating concatenated nucleotide FASTA $fna_unaln..." 
cat "$outdir"/fna_each/*fna > "$fna_unaln"
ls -lh "$fna_unaln"

# Align the sequences
log_time "Running MAFFT to align the nucleotide FASTA..."
runstats mafft \
        --reorder --auto --adjustdirection --leavegappyregion \
        --thread "$threads" --quiet "$fna_unaln" > "$fna_aln"
sed -i -E 's/^>_R_/>/' "$fna_aln" # Fix FASTA header - Remove '_R_' prefixes that MAFFT adds for rev-comp seqs

# Create the tree
log_time "Running IQtree to build a nucleotide tree..."
runstats iqtree \
    -s "$fna_aln" --prefix "$outdir"/trees/"$OG"_nuc \
    $boot_opt -redo -nt "$threads" -ntmax "$threads" -mem "$mem_gb" \
    > "$outdir"/logs/iqtree_"$OG"_nuc.log
log_time "...Done. IQtree log:"
ls -lh "$outdir"/logs/iqtree_"$OG"_nuc.log

# Report
log_time "Done with script orthotree.sh"
log_time "Listing the output trees:"
ls -lh "$nuc_tree" "$prot_tree" "$root_file"
