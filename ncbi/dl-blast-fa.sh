#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-dl-blast-fa-%j.out


# SETUP ------------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "## $0: Download FASTA files from BLAST output."
    echo
    echo "## Syntax: $0 -i <input> -o <output> [-f] [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -i     Input file (Processed BLAST output) (REQUIRED)"
    echo "## -o     Output dir (FASTA) (REQUIRED)"
    echo "## -f     Only download full genomes (default: download all)"
    echo "## Example: $0 -i blast_processed.out -o blast_seqs.fa"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
infile=""                 # BLAST output processed by blast-process.sh
outdir=""                 # Output dir for FASTA files with downloaded seqs
genomes_only=false        # Whether to only download sequences that represent full genomes

## Parse command-line options
while getopts ':i:o:fh' flag; do
    case "${flag}" in
    i) infile="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    f) genomes_only="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script dl-blast-fa.sh..."
date
echo

## Software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/entrez-direct
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08

## Bash strict mode
set -euo pipefail

## Test options
[[ ! -f $infile ]] && echo "## $0: ERROR: Input file $infile does not exist" >&2 && exit 1
[[ $outdir = "" ]] && echo "## $0: ERROR: No output dir (-o) provided" >&2 && exit 1

## Process options
mkdir -p "$outdir"

## Report
echo "## Input file (BLAST output):       $infile"
echo "## Output dir:                      $outdir"
echo "## Download full genomes only?:     $genomes_only"
echo -e "------------------------\n"


# DOWNLOAD FASTA FILES ---------------------------------------------------------
for query in $(cut -f2 "$infile" | sort | uniq); do

    outfile="$outdir"/"$query"_hits.fa

    if [ "$genomes_only" = false ]; then
        accessions=($(awk -v query="$query" '$2 == query' "$infile" | cut -f1 | sort | uniq))
        echo -e "\n## Downloading all ${#accessions[@]} accessions for query $query..."
    else
        accessions=($(awk -v query="$query" '$2 == query' "$infile" | grep "complete genome" | cut -f1 | sort | uniq))
        echo -e "\n## Downloading ${#accessions[@]} full-genome accessions for query $query..."
    fi

    for accession in "${accessions[@]}"; do
        echo "## Accession: $accession" >&2
        esearch -db nuccore -query "$accession" | efetch -format fasta
    done > "$outfile"

done


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Showing output files:"
ls -lh "$outdir"/*_hits.fa

echo -e "\n## Done with script dl-blast-fasta.sh"
date
echo