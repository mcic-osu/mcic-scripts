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
  echo "Syntax: $0 -i <input-dir> -o <output-dir>..."
  echo
  echo "Required options:"
  echo "    -i DIR            Input dir with FASTA files"
  echo "    -o DIR            Output dir"
  echo
  echo "Other options:"
  echo "    -k INTEGER        Kmer size (should be an odd integer)     [default: 31]"
  echo "    -h                Print this help message and exit"
  echo
  echo "Example:              $0 -i data/refgenomes -o results/sourmash -k 29"
  echo "To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo
  echo "Sourmash documentation: https://sourmash.readthedocs.io"
  echo "Sourmash paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6720031/"
  echo
}

## Option defaults
indir=""
outdir=""
kval=31

## Parse command-line options
while getopts ':i:o:k:h' flag; do
  case "${flag}" in
    i) indir="$OPTARG" ;;
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
csv_out="$outdir"/output/refs_init_ani.csv     # ANI matrix in CSV format
cmp_out="$outdir"/output/refs_init_ani.cmp     # ANI matrix in Python format for sourmash plotting
fa_files=( $(find "$indir" -iname '*fasta' -or -iname '*fa' -or -iname '*fna' -or -iname '*fna.gz') )

## Report
echo "## Starting script sourmash_ani.sh"
date
echo
echo "## Dir with input FASTA files:   $indir"
echo "## Output dir:                   $outdir"
echo "## Kmer value:                   $kval"
echo
echo "## FASTA files:"
echo "${fa_files[@]}"
echo -e "--------------------\n"


# RUN SOURMASH -----------------------------------------------------------------
## Create a signature for each FASTA file
echo "## Creating sourmash signatures..."
sourmash sketch dna \
    -p k="$kval" \
    --outdir "$outdir"/signatures \
    "${fa_files[@]}"

## Rename signatures so they have short names
echo -e "--------------------\n"
echo "## Renaming sourmash signatures..."
for sig in "$outdir"/signatures/*sig; do
    newname=$(basename "$sig" .fasta.sig | sed 's/^Spades//')
    newfile="$outdir"/sig_renamed/"$(basename "$sig")"

    sourmash signature rename \
        "$sig" \
        "$newname" \
        -o "$newfile"
done

## Compute ANI values
echo -e "--------------------\n"
echo "## Computing ANI values..."
sourmash compare \
    -k"$kval" \
    --ani \
    --from-file <(ls "$outdir"/sig_renamed/*) \
    --csv "$csv_out" \
    -o "$cmp_out"

## Make ANI plots
echo -e "\n## Creating plots..."
sourmash plot \
    --labels \
    --output-dir "$outdir"/output \
    "$cmp_out"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing the output files"
ls -lh "$outdir"/output
echo -e "\n## Done with script sourmash_ani.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
