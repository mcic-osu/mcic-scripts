#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
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
    echo "  sbatch $0 --entap_dir <EnTap results dir> --asm_1trans <FASTA file> --asm_alltrans <FASTA file> --asm_out <FASTA file> --annot_out <TSV file> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  --entap_dir     <dir>   Dir with EnTap output (input for this script)"
    echo "  --asm_1trans    <file>  Input FASTA file with a single transcript per gene"
    echo "  --asm_alltrans  <file>  Input FASTA file all transcripts per gene"
    echo "  --asm_out       <file>  Output assembly FASTA"
    echo "  --annot_out     <file>  Output (filtered) annotation TSV"
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
Load_seqtk() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqtk
}

Load_bioawk() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate 2>/dev/null; done
    source activate /fs/project/PAS0471/jelmer/conda/bioawk
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
Load_seqtk

# Bash script settings
set -euo pipefail

# Input files
entap_finaldir="$entap_dir"/final_results
entap_final="$entap_finaldir"/final_annotations_no_contam_lvl0.tsv
entap_contam="$entap_finaldir"/final_annotations_contam_lvl0.tsv
frame_dir="$entap_dir"/frame_selection/TransDecoder/processed

# Output files
outdir=$(dirname "$asm_out")
outdir_annot=$(dirname "$annot_out")

not_annot="$outdir"/genes_notannot.txt              # Genes that weren't annotated
gene2trans="$outdir"/gene2trans.tsv                 # Gene-to-transcript map
gene_lengths="$outdir"/gene_lengths.tsv             # Gene lengths for GO
trans_final="$outdir"/transIDs_final.txt            # List of transcripts

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
[[ ! -f "$entap_final" ]] && Die "Input file $entap_final does not exist"
[[ ! -f "$entap_contam" ]] && Die "Input file $entap_contam does not exist"

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
echo
echo "Listing the input files:"
ls -lhd "$entap_dir" 
ls -lh "$asm_1trans" "$asm_alltrans"
echo "=========================================================================="


# ==============================================================================
#                               RUN
# ==============================================================================
# Create the output directory
echo -e "\n# Creating the output directories..."
mkdir -pv "$outdir"/logs "$outdir_annot"


# REMOVE NON-ANNOTATED GENES ---------------------------------------------------
echo -e "\n# Removing non-annotated genes..."

# Make a list of not-annotated genes
join -v 1 \
    <(grep "^>" "$asm_1trans" | sort -u) \
    <(grep "^>" "$entap_finaldir"/final_annotated.faa | sort -u) |
    sed -e 's/>//' -e 's/ .*//' \
    >"$not_annot"

# Exclude not-annotated genes from "final" Entap assembly
grep -v -f "$not_annot" "$entap_final" > "$annot_out"

# Get stats to report
n_in=$(grep -c "^>" "$asm_1trans")
n_notannot=$(wc -l < "$not_annot")
n_out=$(tail -n +2 "$annot_out" | wc -l)


# SUBSET ASSEMBLY WITH ALL ISOFORMS TO KEEP ENTAP SELECTION --------------------
echo -e "\n# Subsetting the assembly with all isoforms..."

# Create a gene-to-transcript map
paste <(grep "^>" "$asm_alltrans" | sed -E 's/>([^ ]+) .*/\1/' | sed -E 's/t[0-9]+//') \
    <(grep "^>" "$asm_alltrans" | sed -E 's/>([^ ]+) .*/\1/') |
    sort -k1,1 >"$gene2trans"

# Create a list of final transcripts to keep
join -t $'\t' "$gene2trans" \
    <(tail -n +2 "$annot_out" | cut -f 1 | sed 's/t1//' | sort) | \
    cut -f 2 >"$trans_final"

# Subset the FASTA file
seqtk subseq "$asm_alltrans" "$trans_final" >"$asm_out"

# Get stats to report
ntrans_out=$(grep -c "^>" "$asm_out")
ntrans_in=$(grep -c "^>" "$asm_alltrans")


# GET GENE LENGTHS FOR GO ------------------------------------------------------
echo -e "\n# Getting the gene lengths..."
Load_bioawk
bioawk -c fastx '{print $name, length($seq)}' "$asm_out" > "$gene_lengths"


# CHECK ENTAP RESULTS ----------------------------------------------------------
if [[ "$check_frames" = true ]]; then
    # Check frame selection
    n_noframe=$(grep "^>" "$frame_dir"/sequences_removed.fnn | sort | uniq | wc -l)
    n_partial=$(grep "^>" "$frame_dir"/partial_genes.faa | sort | uniq | wc -l)
    n_internal=$(grep "^>" "$frame_dir"/internal_genes.faa | sort | uniq | wc -l)
    n_complete=$(grep "^>" "$frame_dir"/complete_genes.faa | sort | uniq | wc -l)
else
    echo "Note: TransDecoder frame selection files are not present."
fi

# Check contamination
n_contam=$(tail -n +2 "$entap_contam" | wc -l)


# REPORT -----------------------------------------------------------------------
echo -e "\n# Done. Summary of the results:"
echo "Number of input genes:                     $n_in"
echo "Number of output genes:                    $n_out"
echo
echo "Number of input transcripts:               $ntrans_in"
echo "Number of output transcripts:              $ntrans_out"
echo
echo "Number of genes marked as contaminant:     $n_contam"
echo "Number of non-annotated genes:             $n_notannot"
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
echo
echo "========================================================================="
echo -e "\n# Listing the output files:"
echo "# Output assembly:"
ls -lh "$asm_out"
echo "# Output annotation:"
ls -lh "$annot_out"
echo "# Output gene-to-transcript map:"
ls -lh "$gene2trans"
echo "# Gene lengths file:"
ls -lh "$gene_lengths"
echo "# List of genes that were not annotated:"
ls -lh "$not_annot"
echo "# List of final transcript IDs:"
ls -lh "$trans_final"

echo -e "\n# Done with script"
date


# SANDBOX ----------------------------------------------------------------------
# Check numbers of genes
# grep -c "^>" "$entap_finaldir"/final_annotated.fnn     # 17,494
# grep -c "^>" "$entap_finaldir"/final_annotations_no_contam.fnn # 36,208
# grep "^>" "$entap_finaldir"/final_annotations_no_contam.fnn | sort | uniq | wc -l # 18,104 -- every gene is present twice in the file...

# Note: Non-contam + contam = every gene, so unannotated must be in those files, too
# tail -n +2 "$entap_final" | wc -l # 18,104
# tail -n +2 "$entap_finaldir"/final_annotations_contam_lvl0.tsv | wc -l    # 3,474 (+18,104=21,578)
