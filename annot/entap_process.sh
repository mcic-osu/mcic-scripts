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
## Help function
Print_help() {
    echo
    echo "======================================================================"
    echo "                            $0"
    echo "                  Process and check EnTAP output"
    echo "======================================================================"
    echo
    echo "USAGE:"
    echo "  sbatch $0 -i <input file> -o <output dir> [...]"
    echo "  bash $0 -h"
    echo
    echo "REQUIRED OPTIONS:"
    echo "  -i/--entap-dir  <dir>   Dir with EnTap output (input for this script)"
    echo "  --in-1trans     <file>  Assembly FASTA file with a single transcript per gene"
    echo "  --in-alltrans   <file>  Assembly FASTA file all transcripts per gene"
    echo "  -o/--outdir    <dir>    Output dir"
    echo "                          Includes an assembly with only transcripts for annotated genes not flagged as contaminants"
    echo
    echo "OTHER KEY OPTIONS:"
    echo "  --more-args     <str>   Quoted string with additional argument(s) to pass to TODO_THIS_SOFTWARE"
    echo
    echo "UTILITY OPTIONS:"
    echo "  --debug                 Run the script in debug mode (print all code)"
    echo "  -h / --help             Print this help message and exit"
    echo
    echo "EXAMPLE COMMANDS:"
    echo "  sbatch $0 -i assembly.fa -c entap_config.ini -d uniprot_sprot.dmnd -o results/entap"
    echo
}

## Load software
Load_seqtk() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/ess/PAS0471/jelmer/conda/seqtk
}

Load_bioawk() {
    module load miniconda3/4.12.0-py39
    [[ -n "$CONDA_SHLVL" ]] && for i in $(seq "${CONDA_SHLVL}"); do source deactivate; done
    source activate /fs/project/PAS0471/jelmer/conda/bioawk
}

## Exit upon error with a message
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
## Option defaults
debug=false


# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
## Placeholder defaults
entap_dir=""
in_1trans=""            # Input assembly with 1 transcript per genes
in_alltrans=""          # Input assembly with all transcripts
outdir=""

## Parse command-line args
all_args="$*"

while [ "$1" != "" ]; do
    case "$1" in
        -i | --entap-dir )  shift && entap_dir=$1 ;;
        --in-1trans )       shift && in_1trans=$1 ;;
        --in-alltrans )     shift && in_alltrans=$1 ;;
        -o | --outdir )     shift && outdir=$1 ;;
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
## In debugging mode, print all commands
[[ "$debug" = true ]] && set -o xtrace

## Load software
Load_seqtk

## Bash script settings
set -euo pipefail

## Input files
entap_findir="$entap_dir"/final_results
entap_fin="$entap_findir"/final_annotations_no_contam_lvl0.tsv
entap_contam="$entap_findir"/final_annotations_contam_lvl0.tsv
frame_dir="$entap_dir"/frame_selection/TransDecoder/processed

## Output files
assembly_out="$outdir"/assembly.fasta
not_annot="$outdir"/genes_notannot.txt                      # Genes that weren't annotated
out_annot="$outdir"/final_annotations_no_contam_lvl0.tsv    # Final annotation after excluding not-annotated genes
gene_trans_map="$outdir"/gene_trans_map.tsv
gene_lengths="$outdir"/gene_lengths.tsv
trans_final="$outdir"/transIDs_final.txt

## Check input
[[ "$entap_dir" = "" ]] && Die "Please specify an input dir with -i/--entap_dir" "$all_args"
[[ "$outdir" = "" ]] && Die "Please specify an output dir with -o/--outdir" "$all_args"
[[ ! -f "$in_1trans" ]] && Die "Input file $in_1trans does not exist"
[[ ! -f "$in_alltrans" ]] && Die "Input file $in_alltrans does not exist"

[[ ! -d "$entap_findir" ]] && Die "Input file $entap_findir does not exist"
[[ ! -f "$entap_fin" ]] && Die "Input file $entap_fin does not exist"
[[ ! -f "$entap_contam" ]] && Die "Input file $entap_contam does not exist"
[[ ! -d "$frame_dir" ]] && Die "Input file $frame_dir does not exist"

## Report
echo
echo "=========================================================================="
echo "                    STARTING SCRIPT TODO_SCRIPTNAME"
date
echo "=========================================================================="
echo "All arguments to this script:             $all_args"
echo "Input dir with EnTAP results:             $entap_dir"
echo "Input assembly - 1 transcript per gene:   $in_1trans"
echo "Input assembly - all transcripts:         $in_alltrans"
echo "Output dir:                               $outdir"
echo "Output assembly:                          $assembly_out"
echo
echo "Listing the input file(s):"
ls -lh "$entap_dir" "$in_1trans" "$in_alltrans"
echo "=========================================================================="


# ==============================================================================
#                               RUN
# ==============================================================================
## Create the output directory
mkdir -p "$outdir"/logs


# REMOVE NON-ANNOTATED GENES ---------------------------------------------------
## Make a list of not-annotated genes
join -v 1 \
    <(grep "^>" "$in_1trans" | sort -u) \
    <(grep "^>" "$entap_findir"/final_annotated.fnn | sort -u) |
    sed -e 's/>//' -e 's/ .*//' \
    >"$not_annot"

## Exclude not-annotated genes from "final" Entap assembly
grep -v -f "$not_annot" "$entap_fin" > "$out_annot"

## Get stats to report
n_in=$(grep -c "^>" "$in_1trans")
n_notannot=$(wc -l < "$not_annot")
n_out=$(tail -n +2 "$out_annot" | wc -l)


# SUBSET ASSEMBLY WITH ALL ISOFORMS TO KEEP ENTAP SELECTION --------------------
## Create a gene-to-transcript map
paste <(grep "^>" "$in_alltrans" | sed -E 's/>([^ ]+) .*/\1/' | sed -E 's/t[0-9]+//') \
    <(grep "^>" "$in_alltrans" | sed -E 's/>([^ ]+) .*/\1/') |
    sort -k1,1 >"$gene_trans_map"

## Create a list of final transcripts to keep
join -t $'\t' "$gene_trans_map" <(tail -n +2 "$out_annot" | cut -f 1 | sed 's/t1//' | sort) | \
    cut -f 2 >"$trans_final"

## Subset the FASTA file
seqtk subseq "$in_alltrans" "$trans_final" >"$assembly_out"

## Get stats to report
ntrans_out=$(grep -c "^>" "$assembly_out")
ntrans_in=$(grep -c "^>" "$in_alltrans")


# GET GENE LENGTHS FOR GO ------------------------------------------------------
Load_bioawk
bioawk -c fastx '{ print $name, length($seq) }' "$assembly_out" > "$gene_lengths"


# CHECK ENTAP RESULTS ----------------------------------------------------------
## Check frame selection
n_noframe=$(grep "^>" "$frame_dir"/sequences_removed.fnn | sort | uniq | wc -l)
n_partial=$(grep "^>" "$frame_dir"/partial_genes.faa | sort | uniq | wc -l)
n_internal=$(grep "^>" "$frame_dir"/internal_genes.faa | sort | uniq | wc -l)
n_complete=$(grep "^>" "$frame_dir"/complete_genes.faa | sort | uniq | wc -l)

## Check contamination
n_contam=$(tail -n +2 "$entap_contam" | wc -l)


# REPORT -----------------------------------------------------------------------
echo "Number of input genes:                     $n_in"
echo "Number of output genes:                    $n_out"
echo
echo "Number of input transcripts:               $ntrans_in"
echo "Number of output transcripts:              $ntrans_out"
echo
echo "Number of genes with no frame:             $n_noframe"
echo "Number of non-annotated genes:             $n_notannot"
echo "Number of genes marked as contaminant:     $n_contam"
echo
echo "Number of partial genes:                   $n_partial"
echo "Number of internal genes:                  $n_internal"
echo "Number of complete genes:                  $n_complete"


# ==============================================================================
#                               WRAP-UP
# ==============================================================================
echo
echo "========================================================================="
echo -e "\n# Listing files in the output dir:"
ls -lhd "$PWD"/"$outdir"/*
echo "# Done with script"
date


# SANDBOX ----------------------------------------------------------------------
## Check numbers of genes
# grep -c "^>" "$entap_findir"/final_annotated.fnn     # 17,494
# grep -c "^>" "$entap_findir"/final_annotations_no_contam.fnn # 36,208
# grep "^>" "$entap_findir"/final_annotations_no_contam.fnn | sort | uniq | wc -l # 18,104 -- every gene is present twice in the file...

## Note: Non-contam + contam = every gene, so unannotated must be in those files, too
# tail -n +2 "$entap_fin" | wc -l # 18,104
# tail -n +2 "$entap_findir"/final_annotations_contam_lvl0.tsv | wc -l    # 3,474 (+18,104=21,578)
