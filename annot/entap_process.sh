#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-entap-process-%j.out

# PARSE OPTIONS ----------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Process and check EnTAP output."
    echo
    echo "Syntax: $0 -d <EnTap-dir> -o <output-dir> -a <assembly-1trans> -A <assembly-alltrans>"
    echo
    echo "Required options:"
    echo "    -i DIR         Dir with EnTap output (input for this script)"
    echo "    -o FILE        Output assembly"
    echo "                   (with only transcripts for annotated genes not flagged as contaminants"
    echo "    -a FILE        Assembly FASTA file with a single transcript per gene"
    echo "    -A FILE        Assembly FASTA file all transcripts per gene"
    echo
    echo "Other options:"
    echo "    -h                Print this help message and exit"
    echo
    echo "Example:    $0 -i assembly.fa -c entap_config.ini -d uniprot_sprot.dmnd -o results/entap"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
entap_dir=""
out_assembly=""
in_1trans=""
in_alltrans=""

## Parse command-line options
while getopts 'i:o:a:A:h' flag; do
    case "${flag}" in
        i) entap_dir="$OPTARG" ;;
        o) out_assembly="$OPTARG" ;;
        a) in_1trans="$OPTARG" ;;
        A) in_alltrans="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/seqtk

## Bash strict mode
set -euo pipefail

## Input files
entap_findir="$entap_dir"/final_results
entap_fin="$entap_findir"/final_annotations_no_contam_lvl0.tsv
entap_contam="$entap_findir"/final_annotations_contam_lvl0.tsv
frame_dir="$entap_dir"/frame_selection/TransDecoder/processed

## Output files
outdir=$(dirname "$out_assembly")
not_annot="$outdir"/genes_notannot.txt                      # Genes that weren't annotated
out_annot="$outdir"/final_annotations_no_contam_lvl0.tsv    # Final annotation after excluding not-annotated genes
gene_trans_map="$outdir"/gene_trans_map.tsv
trans_final="$outdir"/transIDs_final.txt

## Report
echo
echo "## Starting script entap_process.sh"
date
echo
echo "## EnTap dir:                            $entap_dir"
echo "## Output assembly:                      $out_assembly"
echo "## Output annotation:                    $out_annot"
echo -e "------------------------\n"

## Check input
[[ ! -d "$entap_findir" ]] && echo "## ERROR: Dir $entap_findir does not exist" >&2 && exit 1
[[ ! -f "$entap_fin" ]] && echo "## ERROR: File $entap_findir does not exist" >&2 && exit 1
[[ ! -f "$entap_contam" ]] && echo "## ERROR: File $entap_findir does not exist" >&2 && exit 1
[[ ! -d "$frame_dir" ]] && echo "## ERROR: Dir $frame_dir does not exist" >&2 && exit 1

## Create output dir
mkdir -p "$outdir"


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
seqtk subseq "$in_alltrans" "$trans_final" >"$out_assembly"

## Get stats to report
ntrans_out=$(grep -c "^>" "$out_assembly")
ntrans_in=$(grep -c "^>" "$in_alltrans")


# CHECK ENTAP RESULTS ----------------------------------------------------------
## Check frame selection
n_noframe=$(grep "^>" "$frame_dir"/sequences_removed.fnn | sort | uniq | wc -l)
n_partial=$(grep "^>" "$frame_dir"/partial_genes.faa | sort | uniq | wc -l)
n_internal=$(grep "^>" "$frame_dir"/internal_genes.faa | sort | uniq | wc -l)
n_complete=$(grep "^>" "$frame_dir"/complete_genes.faa | sort | uniq | wc -l)

## Check contamination
n_contam=$(tail -n +2 "$entap_contam" | wc -l)


# WRAP-UP ----------------------------------------------------------------------
## Report
echo "## Number of input genes:                            $n_in"
echo "## Number of output genes:                           $n_out"
echo
echo "## Number of input transcripts:                      $ntrans_in"
echo "## Number of output transcripts:                     $ntrans_out"
echo
echo "## Number of genes with no frame:                    $n_noframe"
echo "## Number of non-annotated genes:                    $n_notannot"
echo "## Number of genes marked as contaminant:            $n_contam"
echo
echo "## Number of partial genes:                          $n_partial"
echo "## Number of internal genes:                         $n_internal"
echo "## Number of complete genes:                         $n_complete"

echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script entap_process.sh"
date
echo


# SANDBOX ----------------------------------------------------------------------
## Check numbers of genes
# grep -c "^>" "$entap_findir"/final_annotated.fnn     # 17,494
# grep -c "^>" "$entap_findir"/final_annotations_no_contam.fnn # 36,208
# grep "^>" "$entap_findir"/final_annotations_no_contam.fnn | sort | uniq | wc -l # 18,104 -- every gene is present twice in the file...

## Note: Non-contam + contam = every gene, so unannotated must be in those files, too
# tail -n +2 "$entap_fin" | wc -l # 18,104
# tail -n +2 "$entap_findir"/final_annotations_contam_lvl0.tsv | wc -l    # 3,474 (+18,104=21,578)
