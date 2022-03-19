#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-blast-process-%j.out

# SETUP ------------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "## $0: Process a (fmt-6) BLAST output file."
    echo
    echo "## Syntax: $0 -i <input> -o <output> [ -t <nr-top-hits> ]"
    echo
    echo "## Required options:"
    echo "## -i STRING        Input file (Raw BLAST output)"
    echo "## -o STRING        Output file (Processed BLAST output)"
    echo
    echo "## Other options:"
    echo "## -t INTEGER       Number of top hits [default: 10]"
    echo "## -h               Print this help message and exit"
    echo
    echo "## Example: $0 -q blast.out -o blast_processed.out -t 50"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
blast_out_raw=""        # Raw BLAST output file (= input for this script)
blast_out_proc=""       # Processed BLAST output file (= output of this script)
top_x_hits=10           # Take the top x hits per query (default: 10)

## Parse command-line options
while getopts ':i:o:t:h' flag; do
    case "${flag}" in
    i) blast_out_raw="$OPTARG" ;;
    o) blast_out_proc="$OPTARG" ;;
    t) top_x_hits="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script blast-process.sh..."
date
echo

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/blast-env

## Bash strict mode
#set -euo pipefail

## Process options
outdir=$(dirname "$blast_out_proc")
blast_out_sorted="$outdir"/hits_sorted.txt
blast_out_top="$outdir"/top_"$top_x_hits"_hits.txt
organism_lookup="$outdir"/organism_lookup.txt

## Create output dir, if needed
mkdir -p "$outdir"

## Test input
[[ ! -f $blast_out_raw ]] && echo "## $0: ERROR: Input file $blast_out_raw does not exist" >&2 && exit 1
[[ $blast_out_proc = "" ]] && echo "## $0: ERROR: No output file (-o) provided" >&2 && exit 1

## Report
echo "## BLAST raw output file (input):           $blast_out_raw"
echo "## BLAST processed output file (output):    $blast_out_proc"
echo "## Take the top-x hits, x is:               $top_x_hits"
echo -e "-----------------------\n"


# PROCESS BLAST OUTPUT ---------------------------------------------------------
echo "## Step 1: Sorting BLAST output by goodness of the match"
export LC_ALL=C
sort -k1,1 -k11,11g -k12,12gr -k3,3gr "$blast_out_raw" > "$blast_out_sorted"    # (Sort by: query name (1), then e-value (11), then bitscore (12), then % identical (3))

echo -e "\n## Step 2: Get the top-x hits for each query"
for query in $(cut -f1 "$blast_out_sorted" | sort -u); do
  grep -w -m "$top_x_hits" "$query" "$blast_out_sorted"                         # -w: whole word matching / -m n: stop after n matches 
done > "$blast_out_top"

echo -e "\n## Step 3: Create a lookup table with accession numbers, organism taxonomy, sequence length and full sequence name"
### (See https://www.biostars.org/p/367121/)
### (Accession nr examples for testing: accession=KF772785.1 / accession=MH231153.1)
while read -r line; do

    accession=$(echo "$line" | cut -f 2)
    esearch_result=$(esearch -db nuccore -query "$accession" </dev/null)
    
    ### Get the taxon name
    taxon=$(echo "$esearch_result" |
                elink -target taxonomy |
                efetch -format native -mode xml |
                grep "ScientificName" | head -n1 |
                sed -E 's@</?ScientificName>@@g' | sed -e 's/^[ \t]*//')
    
    ### Get the sequence length
    seq_len=$(echo "$esearch_result" |
                efetch -format native -mode xml |
                grep "Seq-inst_length" | head -n1 |
                sed -E 's@</?Seq-inst_length>@@g' | sed -e 's/^[ \t]*//')

    ### Get the full name of the sequence / NCBI entry
    seq_title=$(echo "$esearch_result" |
                    efetch -format native -mode xml |
                    grep "Seqdesc_title" | head -n1 |
                    sed -E 's@</?Seqdesc_title>@@g' | sed -e 's/^[ \t]*//')

    echo -e "${accession}\t${taxon}\t${seq_len}\t${seq_title}"
    echo -e "${accession}\t${taxon}\t${seq_len}\t${seq_title}" >&2

done < "$blast_out_top" > "$organism_lookup"

echo -e "\n## Step 4: Merge the BLAST output table with the taxonomy lookup table"
join -t $'\t' -1 2 -2 1 "$blast_out_top" "$organism_lookup" > "$blast_out_proc"


# WRAP UP ----------------------------------------------------------------------
echo -e "\n\n## Listing the output file:"
ls -lh "$blast_out_proc"

echo -e "\n## Done with script blast-process.sh"
date
echo
