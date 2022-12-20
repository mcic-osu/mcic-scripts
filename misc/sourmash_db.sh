#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=sourmash_db
#SBATCH --output=slurm-sourmash-db-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
    echo
    echo "$0: Create a custom sourmash database."
    echo
    echo "Syntax: $0 -i <input-dir> -o <output-dir> ..."
    echo
    echo "Required options:"
    echo "    -i DIR            Input dir with FASTA files (extensions .fa/.fasta/.fna)"
    echo "    -o DIR            Output dir"
    echo
    echo "Other options:"
    echo "    -d STRING         Database name"
    echo "    -k INTEGER        Kmer size (should be an odd integer)     [default: 31]"
    echo "    -h                Print this help message and exit"
    echo
    echo "Example:              $0 -i data/refgenomes -o results/sourmash/db -d mydb -k 29"
    echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
    echo "Sourmash documentation: https://sourmash.readthedocs.io/en/latest/tutorial-basic.html#make-and-search-a-database-quickly"
    echo "Sourmash paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/"
    echo
}

## Option defaults
indir=""
outdir=""
kval=31
db_name=smash_db

## Parse command-line options
while getopts ':i:o:d:k:h' flag; do
    case "${flag}" in
        i) indir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        d) db_name="$OPTARG" ;;
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
[[ ! $indir =~ ^/ ]] && indir="$PWD"/"$indir"

## Report
echo "## Starting script sourmash_db.sh"
date
echo
echo "## Input dir:                    $indir"
echo "## Output dir:                   $outdir"
echo "## Database name:                $db_name"
echo "## Kmer value:                   $kval"
echo -e "--------------------\n"


# RUN SOURMASH -----------------------------------------------------------------
## Move to output dir
cd "$outdir" || exit

## Create signatures (sketches) for each input file
for fa in $(find "$indir" -iname '*fasta' -or -iname '*fa' -or -iname '*fna' -or -iname '*fna.gz'); do
    if [[ ! -f $(basename "$fa").sig ]]; then
        sourmash sketch dna -p k="$kval" "$fa"
    else
        echo "Signature file already exists for $fa"
    fi
done

## Create the database
sourmash index -k"$kval" "$db_name" ./*.sig


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh
echo -e "\n## Done with script sourmash_db.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
