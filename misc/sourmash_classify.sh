#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sourmash_classify
#SBATCH --output=slurm-sourmash-classify-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Search for a query in a sourmash database using LCA classification."
    echo
    echo "Syntax: $0 -i <input-dir> -o <output-dir> ..."
    echo
    echo "Required options:"
    echo "    -i FILE           Input FASTA file"
    echo "    -o DIR            Output dir"
    echo
    echo "Other options:"
    echo "    -d FILE           Path to a .lca.json.gz sourmash database"
    echo "                      [default: download GTDB database]"
    echo "    -D DIR            Directory to download GTDB database to (use either -d or -D)"
    echo "                      [default: download to output dir]"
    echo "    -k INTEGER        Kmer size (21, 31, or 51)                [default: 31]"
    echo "    -h                Print this help message and exit"
    echo
    echo "Example:              $0 -i data/refgenomes -o results/sourmash/db -d mydb -k 29"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "Sourmash documentation: https://sourmash.readthedocs.io"
    echo "Sourmash paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/"
    echo
}

## Option defaults
fa_in=""
db=""
db_dir=""
outdir=""
kval=31

## Parse command-line options
while getopts ':i:o:d:D:k:h' flag; do
    case "${flag}" in
        i) fa_in="$OPTARG" ;;
        d) db="$OPTARG" ;;
        D) db_dir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        k) kval="$OPTARG" ;;
        h) Help && exit 0 ;;
        \?) echo -e "\n## $0: ERROR: Invalid option -$OPTARG\n\n" >&2 && exit 1 ;;
        :) echo -e "\n## $0: ERROR: Option -$OPTARG requires an argument\n\n" >&2 && exit 1 ;;
    esac
done


# SETUP ------------------------------------------------------------------------
## Software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/sourmash

## Bash strict mode
set -euo pipefail

## Output
mkdir -p "$outdir"

## If needed, make dirs absolute because we have to move into the outdir
[[ ! $fa_in =~ ^/ ]] && fa_in="$PWD"/"$fa_in"
[[ ! $db =~ ^/ ]] && db="$PWD"/"$db"
[[ $db_dir != "" ]] && [[ ! $db_dir =~ ^/ ]] && db_dir="$PWD"/"$db_dir"
[[ ! $outdir =~ ^/ ]] && outdir="$PWD"/"$outdir"

## Process args
[[ "$db_dir" = "" && "$db" = "" ]] && db_dir="$outdir"
sig_in=$(basename "$fa_in").sig
smpID=$(basename "$fa_in" | sed -E 's/\.fn?as?t?a?//')

## Report
echo "## Starting script sourmash_classify.sh"
date
echo
echo "## Input FASTA file:             $fa_in"
[[ "$db" != "" ]] && echo "## Database:                     $db"
[[ "$db_dir" != "" ]] && echo "## Database dir:                 $db_dir"
echo "## Output dir:                   $outdir"
echo "## Kmer value:                   $kval"
echo "## Sample ID:                    $smpID"
echo -e "--------------------\n"


# RUN SOURMASH -----------------------------------------------------------------
## Move to output dir
cd "$outdir" || exit

## Create a signature for the FASTA file
if [[ ! -f "$sig_in" ]]; then
    echo -e "## Create sourmash signature for query file... ($sig_in)"
    sourmash sketch dna -p k="$kval" "$fa_in"
    echo
else
    echo -e "## Sourmash signature file for query already exists ($sig_in)\n"
fi

## Download the database - https://sourmash.readthedocs.io/en/latest/databases.html
if [[ "$db_dir" != "" ]]; then
    
    db="$db_dir"/gtdb-rs207.genomic.k"$kval".lca.json.gz
    
    if [[ ! -f "$db" ]]; then
        echo "## Downloading database for k=$kval"
        echo "## (See https://sourmash.readthedocs.io/en/latest/databases.html#gtdb-all-genomes-258k)"
        [[ "$kval" = 21 ]] && curl -JL -o "$db" https://osf.io/hm3c4/download
        [[ "$kval" = 31 ]] && curl -JL -o "$db" https://osf.io/tf3ah/download
        [[ "$kval" = 51 ]] && curl -JL -o "$db" https://osf.io/3cdp6/download
        [[ ! -f "$db" ]] && echo "ERROR: Downloaded DB file does not exist/have expected name $db" >&2 & exit 1
        echo
    fi
fi

## Run the search
echo -e "## Now running sourmash classify..."
sourmash lca classify \
    --db "$db" \
    --query "$sig_in" |
    tee >"$smpID".txt


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Showing the main output file:"
cat "$smpID".txt
echo -e "\n## Done with script sourmash_classify.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
