#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=entap_process
#SBATCH --output=slurm-entap_process-%j.out

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
# Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  Process and check EnTAP output"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 --entap_dir <EnTap results dir> --asm_1trans <FASTA file> \\"
    echo "            --asm_alltrans <FASTA file> --asm_out <FASTA file> --annot_out <TSV file> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --entap_dir     <dir>   Dir with EnTap output (input for this script)"
    echo "  --asm_1trans    <file>  Input nucleotide FASTA assembly with a single transcript per gene"
    echo "  --asm_alltrans  <file>  Input nucleotide FASTA assembly with all transcripts"
    echo "  --asm_out       <file>  Output nucleotide assembly FASTA (will have all transcripts)"
    echo "  --annot_out     <file>  Output (filtered) annotation TSV"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --eggnog_contam <str>   Quoted, comma-separated list of EggNOG taxa considered contaminants"
    echo "                            - Default: 'Animals,Fungi,Bacteria,Arthropoda,Opisthokonts,Mammals,Fishes,Aves,Archaea,Nematodes'"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h / --help             Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 --entap_dir results/entap --asm_1trans results/assembly.fa --asm_alltrans results/assembly_alltrans.fa \ "
    echo "     --asm_out results/entap_proc/assembly.fa --annot_out results/entap_proc/annot.tsv"
    echo
}

# Load software
Load_software() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqkit
}

# Exit upon error with a message
Die() {
    error_message=${1}
    error_args=${2-none}
    
    echo >&2
    echo "=====================================================================" >&2
    printf "$0: ERROR: %s\n" "$error_message" >&2
    echo -e "\nFor help, run this script with the '-h' option" >&2
    echo "For example, 'bash mcic-scripts/qc/fastqc.sh -h'" >&2
    if [[ "$error_args" != "none" ]]; then
        echo -e "\nArguments passed to the script:" >&2
        echo "$error_args" >&2
    fi
    echo -e "\nEXITING..." >&2
    echo "=====================================================================" >&2
    echo >&2
    exit 1
}


# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Option defaults
eggnog_contam="Animals,Fungi,Bacteria,Arthropoda,Opisthokonts,Mammals,Fishes,Aves,Archaea,Nematodes"

check_frames=true
debug=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Placeholder defaults
entap_dir=""
asm_1trans=""            # Input assembly with 1 transcript per genes
asm_alltrans=""          # Input assembly with all transcripts
asm_out=""
annot_out=""

# Parse command-line args
all_args="$*"
while [ "$1" != "" ]; do
    case "$1" in
        --entap_dir )       shift && entap_dir=$1 ;;
        --asm_1trans )      shift && asm_1trans=$1 ;;
        --asm_alltrans )    shift && asm_alltrans=$1 ;;
        --asm_out )         shift && asm_out=$1 ;;
        --annot_out )       shift && annot_out=$1 ;;
        --eggnog_contam )   shift && eggnog_contam=$1 ;;
        -v | --version )    Print_version; exit 0 ;;
        -h )                Print_help; exit 0 ;;
        --help )            Print_help_program; exit 0;;
        --debug )           debug=true ;;
        * )                 Die "Invalid option $1" "$all_args" ;;
    esac
    shift
done


# ==============================================================================
#                          OTHER SETUP
# ==============================================================================
# In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

# Load software
Load_software

# Bash script settings
set -euo pipefail

# Input files
entap_finaldir="$entap_dir"/final_results
annot_in="$entap_finaldir"/final_annotations_no_contam_lvl0.tsv
simsearch_contam="$entap_finaldir"/final_annotations_contam_lvl0.tsv
frame_dir="$entap_dir"/frame_selection/TransDecoder/processed

# Output dirs
outdir=$(dirname "$asm_out")
outdir_annot=$(dirname "$annot_out")

# Output files
not_annot_ids="$outdir"/notannot_ids.txt                # Genes that weren't annotated
gene2trans="$outdir"/gene2trans.tsv                 # Gene-to-transcript map
gene_lengths="$outdir"/gene_lengths.tsv             # Gene lengths for GO
trans_final="$outdir"/transIDs_final.txt            # List of transcripts
simsearch_contam_ids="$outdir"/simsearch_contam_ids.txt # IDs of transcripts IDed as contaminant by sim. search
eggnog_contam_ids="$outdir"/eggnog_contam_ids.txt   # IDs of transcripts IDed as contaminant by EggNOG
eggnog_contam_taxa="$outdir"/eggnog_contam_taxa.txt # List of taxa considered contaminants

asm_aa_in="$outdir"/assembly_1trans_in.faa
asm_aa_out="$outdir"/assembly_1trans.faa
asm_aa_notannot="$outdir"/assembly_notannot.faa
asm_nuc_notannot="$outdir"/assembly_notannot.fna

asm_nuc_contam_simsearch="$outdir"/assembly_simsearchcontam.fna
asm_aa_contam_simsearch="$outdir"/assembly_simsearchcontam.faa
asm_nuc_contam_eggnog="$outdir"/assembly_eggnogcontam.fna
asm_aa_contam_eggnog="$outdir"/assembly_eggnogcontam.faa

simsearchNo_eggnogYes_ids="$outdir"/nosimsearch_buteggnog_ids.txt
no_simsearch_ids="$outdir"/nosimsearch_ids.txt
asm_aa_nosimsearch="$outdir"/assembly_nosimsearch.faa
asm_nuc_nosimsearch="$outdir"/assembly_nosimsearch.fna

# Output files - temporary
annot_intermed="$outdir"/annot_intermed.tsv         # Before eggnog-contam. removal
annot_out_tmp="$annot_out".tmp

# Check if TransDecoder output is there
[[ ! -d "$frame_dir" ]] && check_frames=false

# Check input
[[ "$entap_dir" = "" ]] && Die "Please specify an dir with EnTap results with --entap_dir" "$all_args"
[[ "$asm_1trans" = "" ]] && Die "Please specify an input assembly with --asm_1trans" "$all_args"
[[ "$asm_alltrans" = "" ]] && Die "Please specify an input assembly with --asm_alltrans" "$all_args"
[[ "$asm_out" = "" ]] && Die "Please specify an output assembly --asm_out" "$all_args"
[[ "$annot_out" = "" ]] && Die "Please specify an output annotation with --annot_out" "$all_args"

[[ ! -d "$entap_dir" ]] && Die "Input dir $entap_dir does not exist"
[[ ! -d "$entap_finaldir" ]] && Die "Input dir $entap_finaldir does not exist"
[[ ! -f "$asm_1trans" ]] && Die "Input file $asm_1trans does not exist"
[[ ! -f "$asm_alltrans" ]] && Die "Input file $asm_alltrans does not exist"
[[ ! -f "$annot_in" ]] && Die "Input file $annot_in does not exist"
[[ ! -f "$simsearch_contam" ]] && Die "Input file $simsearch_contam does not exist"

# Get stats to report
n_genes_in=$(grep -c "^>" "$asm_1trans")
n_trans_in=$(grep -c "^>" "$asm_alltrans")

# Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT ENTAP_PROCESS.SH"
date
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input dir with EnTAP results:             $entap_dir"
echo "Input assembly - 1 transcript per gene:   $asm_1trans"
echo "Input assembly - all transcripts:         $asm_alltrans"
echo "Output assembly:                          $asm_out"
echo "Output annotation:                        $annot_out"
echo
echo "Nr genes in the input assembly:           $n_genes_in"
echo "Nr transcripts in the input assembly:     $n_trans_in"
echo
echo "# Listing the input files:"
ls -lhd "$entap_dir" 
ls -lh "$asm_1trans" "$asm_alltrans"
echo "=========================================================================="


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
mkdir -pv "$outdir"/logs "$outdir_annot"

# Get the annotation header
head -n 1 "$annot_in" > "$outdir"/annot_header.txt


# DEDUPLICATE ENTAP'S AA FASTAS ------------------------------------------------
echo -e "\n===================================================================="
echo "# Removing duplicates from EnTAP's AA FASTA files..."

# Contaminated fasta
seqkit rmdup -n "$entap_finaldir"/final_annotations_contam.faa > "$asm_aa_contam_simsearch"

# Input aa FASTA 
seqkit rmdup -n "$entap_finaldir"/final_annotations.faa > "$asm_aa_in"

# Unannotated aa FASTA (Probably does not contain duplicates anyway)
seqkit rmdup -n "$entap_finaldir"/final_unannotated.faa > "$asm_aa_notannot"

echo "# Listing the output files:"
ls -lh "$asm_aa_contam_simsearch" "$asm_aa_in" "$asm_aa_notannot"


# CHECK SIM-SEARCH CONTAMINANTS ------------------------------------------------
echo -e "\n===================================================================="
echo "# Checking similarity search contaminants..."

# Check how many transcripts were identified as contaminants by similarity search
n_contam=$(tail -n +2 "$simsearch_contam" | wc -l)
echo "There are $n_contam contaminants listed in the file $simsearch_contam" 

# Create a file with just the IDs of contaminant transcripts
tail -n +2 "$simsearch_contam" | cut -f 1 > "$simsearch_contam_ids"

# Create a assemblies with only contaminants
seqkit grep -f "$simsearch_contam_ids" "$asm_1trans" > "$asm_nuc_contam_simsearch"
seqkit grep -f "$simsearch_contam_ids" "$asm_aa_in" > "$asm_aa_contam_simsearch"

# Report
n_nuc=$(grep -c "^>" < "$asm_nuc_contam_simsearch")
n_aa=$(grep -c "^>" < "$asm_aa_contam_simsearch")
echo -e "\n# Assemblies with only contaminants ($n_nuc / $n_aa genes):"
ls -lh "$asm_nuc_contam_simsearch" "$asm_aa_contam_simsearch"


# PROCESS NON-ANNOTATED GENES ---------------------------------------------------
echo -e "\n===================================================================="
echo "# Processing non-annotated genes..."

# A) Make a list of not-annotated genes
grep "^>" "$entap_finaldir"/final_unannotated.faa | sed 's/>//' | sort > "$not_annot_ids"

n_notannot=$(wc -l < "$not_annot_ids")
echo -e "\n# Not-annotated gene-list file ($n_notannot genes):"
ls -lh "$not_annot_ids"
echo -e "\nNumber of non-annotated genes: $n_notannot"

# B) Exclude not-annotated genes from "final" Entap assembly
head -n 1 "$annot_in" > "$annot_intermed"
join -t $'\t' -v 2 "$not_annot_ids" <(tail -n+2 "$annot_in" | sort) >> "$annot_intermed"

n_intermed=$(wc -l < "$annot_intermed")
echo -e "\n# Intermediate annotation file ($n_intermed genes):"
ls -lh "$annot_intermed"

# C) Get a nucleotide assembly with only non-annotated genes
seqkit grep -f "$not_annot_ids" "$asm_1trans" > "$asm_nuc_notannot"

n_notannot2=$(grep -c "^>" "$asm_nuc_notannot")
echo -e "\n# Nucleotide assembly with not-annotated transcripts ($n_notannot2 genes):"
ls -lh "$asm_nuc_notannot"

# D) Get AA and nuc assemblies with genes with no-*simsearch*-annotation only, for re-annotation
tail -n +2 "$annot_intermed" | awk -F"\t" '$2 == ""' | cut -f 1 > "$simsearchNo_eggnogYes_ids"
cat "$not_annot_ids" "$simsearchNo_eggnogYes_ids" | sort -u > "$no_simsearch_ids"
#(Add 'frame' because for Entapnf annotation, the header line needs to have a space and a second word:)
seqkit grep -f "$no_simsearch_ids" "$asm_aa_in" | sed '/^>/s/$/ frame/' > "$asm_aa_nosimsearch"
seqkit grep -f "$no_simsearch_ids" "$asm_1trans" > "$asm_nuc_nosimsearch"

n_nosim=$(grep -c ">" "$asm_aa_nosimsearch")
echo -e "\n# AA assembly with no-simsearch-annotation ($n_nosim genes)"
echo "  (This assembly contains: A. genes w/ no annot. at all; B. genes w/ eggnog but no sim.search annot.)"
ls -lh "$asm_aa_nosimsearch"


# REMOVE EGGNOG-IDENTIFIED CONTAMINANTS ----------------------------------------
echo -e "\n===================================================================="
echo "# Removing genes with no DIAMOND annotation and EggNOG annotation from contaminant taxa..."

# Write Eggnog contaminant taxa to file
echo "$eggnog_contam" | tr "," "\n" > "$eggnog_contam_taxa"

echo -e "\n# Eggnog contaminant taxa:"
cat -n "$eggnog_contam_taxa"

# A) Get the column number which has the EggNOG tax. annotations
eggnog_tax_col=$(head -n 1 "$annot_intermed" | tr "\t" "\n" | cat -n | grep "EggNOG Tax Scope$" | awk '{print $1}')

# B) Identify contaminants
tail -n +2 "$annot_intermed" |
    awk -F"\t" '$2 == ""' |
    awk -F"\t" -v fcol="$eggnog_tax_col" '{print $1 "\t" $fcol}' |
    awk 'NF>1' |
    grep -f "$eggnog_contam_taxa" |
    cut -d " " -f 1 |
    sort \
    > "$eggnog_contam_ids"

n_contam_eggnog=$(wc -l < "$eggnog_contam_ids")
echo -e "\n# Eggnog contaminant IDs file ($n_contam_eggnog genes)"
ls -lh "$eggnog_contam_ids"

# C) Remove contaminants
join -t $'\t' -v 2 "$eggnog_contam_ids" <(tail -n +2 "$annot_intermed") > "$annot_out_tmp"

# D) Create FASTAs with contaminants
seqkit grep -f <(cut -f 1 "$eggnog_contam_ids") "$asm_aa_in" > "$asm_aa_contam_eggnog"
seqkit grep -f <(cut -f 1 "$eggnog_contam_ids") "$asm_1trans" > "$asm_nuc_contam_eggnog"

# Report
n_aa=$(grep -c "^>" < "$asm_aa_contam_eggnog")
n_nuc=$(grep -c "^>" < "$asm_nuc_contam_eggnog")
n_out=$(tail -n +2 "$annot_out_tmp" | wc -l)
echo -e "\n# Output annotation (temporary) ($n_out genes):"
ls -lh "$annot_out_tmp"
echo -e "\n# Assemblies with only contaminants ($n_aa / $n_nuc genes):"
ls -lh "$asm_aa_contam_eggnog" "$asm_nuc_contam_eggnog"


# SUBSET ASSEMBLY WITH ALL ISOFORMS TO KEEP ENTAP SELECTION --------------------
echo -e "\n===================================================================="
echo "# Subsetting the assembly with all isoforms..."

# Create a gene-to-transcript map
paste <(grep "^>" "$asm_alltrans" | awk '{print $1}' | sed 's/>//' | sed -E 's/t[0-9]+$//') \
    <(grep "^>" "$asm_alltrans" | awk '{print $1}' | sed 's/>//') |
    sort -k1,1 >"$gene2trans"

echo -e "\n# Gene-to-transcript map:"
ls -lh "$gene2trans"

# Create a list of final transcripts to keep
join -t $'\t' "$gene2trans" \
    <(tail -n +2 "$annot_out_tmp" | cut -f 1 | sed 's/t1//' | sort) | \
    cut -f 2 >"$trans_final"

ntrans_list=$(wc -l < "$trans_final")
echo -e "\n# Final list of transcripts ($ntrans_list transcripts):"
ls -lh "$trans_final"

# Subset the nucleotide FASTA file
seqkit grep -f "$trans_final" "$asm_alltrans" > "$asm_out"

# Subset the aa FASTA file
seqkit grep -f "$trans_final" "$asm_aa_in" > "$asm_aa_out"

ntrans_nuc_out=$(grep -c "^>" "$asm_out")
ntrans_aa_out=$(grep -c "^>" "$asm_aa_out")
echo -e "\n# Output nuc assembly with all transcripts ($ntrans_nuc_out transcripts):"
ls -lh $"$asm_out"
echo -e "\n# Output AA assembly with 1 transcript per gene ($ntrans_aa_out transcripts):"
ls -lh "$asm_aa_out"


# FIXING THE OUTPUT ANNOTATION FILE --------------------------------------------
echo -e "\n===================================================================="
echo "# Finalizing the output annotation file..."

# Get the header
cat "$outdir"/annot_header.txt > "$annot_out"

# Remove the full paths from the DIAMOND DB origin column
awk -F"\t" -v OFS="\t" '{sub(/.*\//, "", $16)}1' "$annot_out_tmp" >> "$annot_out" 

echo -e "\n# Final annotation file:"
ls -lh "$annot_out"  


# COUNTING ANNOTATION SOURCES --------------------------------------------------
echo -e "\n===================================================================="
echo "# Counting annotation sources..."
tail -n +2 "$annot_out" | awk -F"\t" '{print $16}' | sort | uniq -c

# Eggnog Tax
echo -e "\n# Counting EggNOG tax. assignments:"
tail -n +2 "$annot_out" | awk -F"\t" -v fcol="$eggnog_tax_col" '{print $fcol}' | sort | uniq -c

# Pathways
eggnog_kegg_col=$(head -n 1 "$annot_out" | tr "\t" "\n" | cat -n | grep "EggNOG KEGG Terms" | awk '{print $1}')
n_eggnog_kegg=$(tail -n +2 "$annot_out" | awk -F"\t" -v fcol="$eggnog_kegg_col" '$fcol != ""' | wc -l)
n_eggnog_go=$(tail -n +2 "$annot_out" | awk -F"\t" -v fcol="$eggnog_kegg_col" '$(fcol+1) != "" || $(fcol+2) != "" || $(fcol+3) != ""' | wc -l)


# GET GENE LENGTHS FOR GO ------------------------------------------------------
echo -e "\n===================================================================="
echo "# Getting gene lengths..."
seqkit fx2tab --length --name "$asm_out" > "$gene_lengths"

echo -e "\n# Gene length file:"
ls -lh $"$gene_lengths"


# CHECK ENTAP RESULTS ----------------------------------------------------------
echo -e "\n===================================================================="
if [[ "$check_frames" = true ]]; then
    # Check frame selection
    n_noframe=$(grep "^>" "$frame_dir"/sequences_removed.fnn | sort | uniq | wc -l)
    n_partial=$(grep "^>" "$frame_dir"/partial_genes.faa | sort | uniq | wc -l)
    n_internal=$(grep "^>" "$frame_dir"/internal_genes.faa | sort | uniq | wc -l)
    n_complete=$(grep "^>" "$frame_dir"/complete_genes.faa | sort | uniq | wc -l)
else
    echo "Note: TransDecoder frame selection files are not present."
fi


# REMOVE UNNECESSARY FILES -----------------------------------------------------
echo -e "\n===================================================================="
echo "# Removing unnecessary files..."
rm -v "$annot_out_tmp" "$annot_intermed" "$outdir"/annot_header.txt


# REPORT -----------------------------------------------------------------------
echo -e "\n===================================================================="
echo "Summary of the results:"
echo "Number of input genes:                     $n_genes_in"
echo "Number of output genes:                    $n_out"
echo
echo "Number of input transcripts:               $n_trans_in"
echo "Number of output transcripts:              $ntrans_nuc_out"
echo
echo "Number of contaminant genes - sim.search:  $n_contam"
echo "Number of contaminant genes - eggnog:      $n_contam_eggnog"
echo "Number of non-annotated genes:             $n_notannot"
echo
echo "Nr genes with EggNOG GO terms:             $n_eggnog_go"
echo "Nr genes with EggNOG KEGG terms:           $n_eggnog_kegg"
if [[ "$check_frames" = true ]]; then
    echo
    echo "Number of genes with no frame:             $n_noframe"
    echo "Number of partial genes:                   $n_partial"
    echo "Number of internal genes:                  $n_internal"
    echo "Number of complete genes:                  $n_complete"
fi


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo -e "\n===================================================================="
echo -e "\n# Listing the output files:"
ls -lhd "$PWD"/"$outdir"/*
echo -e "\n# Done with script"
date


# SANDBOX ----------------------------------------------------------------------
# Check numbers of genes
# grep -c "^>" "$entap_finaldir"/final_annotated.fnn     # 17,494
# grep -c "^>" "$entap_finaldir"/final_annotations_no_contam.fnn # 36,208
# grep "^>" "$entap_finaldir"/final_annotations_no_contam.fnn | sort | uniq | wc -l # 18,104 -- every gene is present twice in the file...

# Note: Non-contam + contam = every gene, so unannotated must be in those files, too
# tail -n +2 "$annot_in" | wc -l # 18,104
# tail -n +2 "$entap_finaldir"/final_annotations_contam_lvl0.tsv | wc -l    # 3,474 (+18,104=21,578)
