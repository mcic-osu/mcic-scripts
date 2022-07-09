#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sourmash_ani
#SBATCH --output=slurm-sourmash-ani-%j.out

# PARSE ARGUMENTS --------------------------------------------------------------
## Help function
Help() {
  echo
  echo "$0: Run an ANI (Average Nucleotide Identity) analysis using sourmash."
  echo
  echo "Syntax: $0 -i <input-dir> -o <output-dir> ..."
  echo
  echo "Required options:"
  echo "    -i FILE           Input FASTA file"
  echo "    -d FILE           Path to a sourmash database"
  echo "    -o DIR            Output dir"
  echo
  echo "Other options:"
  echo "    -k INTEGER        Kmer size (should be an odd integer)     [default: 31]"
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
outdir=""
kval=31

## Parse command-line options
while getopts ':i:o:d:k:h' flag; do
  case "${flag}" in
    i) fa_in="$OPTARG" ;;
    d) db="$OPTARG" ;;
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
source activate /fs/project/PAS0471/jelmer/conda/sourmash-4.4.0

## Bash strict mode
set -euo pipefail

## Create output dirs
mkdir -p "$outdir"/signatures "$outdir"/sig_renamed "$outdir"/output

## Process args
sig_in=$(basename "$fa_in").sig
smpID=$(basename "$fa_in" | sed -E 's/\.fn?as?t?a?//')

## Report
echo "## Starting script sourmash_ani.sh"
date
echo
echo "## Input FASTA file:             $fa_in"
echo "## Database:                     $db"
echo "## Output dir:                   $outdir"
echo "## Kmer value:                   $kval"
echo "## Sample ID:                    $smpID"
echo -e "--------------------\n"


# RUN SOURMASH -----------------------------------------------------------------
## Create a signature for each FASTA file
if [[ ! -f "$sig_in" ]]; then
    echo -e "\n## Create sourmash signature for query file..."
    sourmash sketch dna -p k="$kval" "$fa_in"
else
    echo -e "\n## Sourmash signature for query file already exists"
fi

## Rename signatures so they have short names
for sig in "$outdir"/signatures/*sig; do
    newname=$(basename "$sig" .fasta.sig)
    sourmash signature rename "$sig" "$newname" -o "$outdir"/sig_renamed/"$(basename "$sig")"
done

## Compute ANI values
sourmash compare -k31 --ani --from-file <(ls "$outdir"/sig_renamed/*) \
    --csv "$outdir"/output/refs_init_ani.csv -o "$outdir"/output/refs_init_ani.cmp

## Make ANI plots
sourmash plot --labels --output-dir "$outdir"/output "$outdir"/output/refs_init_ani.cmp


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing the output files"
ls -lh "$outdir"/output
echo -e "\n## Done with script sourmash_ani.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
