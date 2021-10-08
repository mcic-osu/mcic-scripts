#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-trimmomatic-%j.out

# Software
module load python/3.6-conda5.2
source activate trimgalore-env

# Bash strict settings
set -euo pipefail

## Help
Help() {
  # Display Help
  echo
  echo "## $0: Run TrimGalore for a FASTQ file."
  echo
  echo "## Syntax: $0 -i <FASTQ input file> -o <FASTQ output dir> -O <FastQC output dir> -q <min seq qual>  -q <min seq len> [-h]"
  echo "## Options:"
  echo "## -h     Print help."
  echo "## -i     Input file (REQUIRED)"
  echo "## -o     Output dir (REQUIRED)"
  echo "## -f     Forward primer (REQUIRED)"
  echo "## -r     Reverse primer (REQUIRED)"
  echo "## -d     Don't discard untrimmed (default: discard)"
  echo "## Example: $0 -i data/fastq/raw -o data/fastq/trimmed -f GAGTGYCAGCMGCCGCGGTAA -r ACGGACTACNVGGGTWTCTAAT [-h]"
  echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
  echo "## The script will compute and use the reverse complements of both primers."
  echo
}

## Option defaults
qual=0    # = no quality trimming
len=20    # = TrimGalore default

# Get command-line options:
while getopts ':i:o:O:q:l:h' flag; do
  case "${flag}" in
  i) fastq_in="$OPTARG" ;;
  o) outdir_trim="$OPTARG" ;;
  O) outdir_fastqc="$OPTARG" ;;
  q) qual="$OPTARG" ;;
  l) len="$OPTARG" ;;
  h) Help && exit 0 ;;
  \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
  :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Process variables and args
mkdir -p "$outdir_trim"
mkdir -p "$outdir_fastqc"

n_threads="$SLURM_CPUS_ON_NODE"

## Report
echo "## Starting script trimgalore.sh"
date
echo "## Input file:                    $fastq_in"
echo "## Output dir - trimmed FASTQs:   $outdir_trim"
echo "## Output dir - FastQC:           $outdir_fastqc"
echo "## Sequence quality threshold:    $qual"
echo "## Minimum sequence length:       $len"
echo -e "---------------------------\n\n"

## Run Trim-Galore
trim_galore \
    --quality "$qual" --length "$len" \
    --gzip -j "$n_threads" \
    --output_dir "$outdir_trim" \
    --fastqc --fastqc_args "-t $n_threads --outdir $outdir_fastqc" \
    "$fastq_in"

## Move FASTQ files
file_id=$(basename "$fastq_in" .fastq.gz)
mv "$outdir_trim"/"$file_id"_trimmed.fq.gz "$outdir_trim"/"$file_id"_trimmed.fastq.gz

## Report
echo -e "\nDone with script trimgalore.sh"
date