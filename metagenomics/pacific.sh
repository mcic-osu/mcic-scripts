#!/bin/bash

#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-pacific-%j.out

# Software:
PACIFIC_DIR=software/pacific
PACIFIC="$PACIFIC_DIR"/PACIFIC.py

source ~/.bashrc
source activate pacific-env

# Strict settings:
set -euo pipefail

## Command-line args:
infile=$1
outdir_base=$2
refseq_dir=$3

file_ID=$(basename "$infile" .fastq.gz)
outdir_full="$outdir_base"/"$file_ID"

mkdir -p "$outdir_full"

## Report:
echo "## Starting script pacific.sh..."
date
echo "## Input FASTQ file: $infile"
echo "## Refseq dir: $refseq_dir"
echo "## Output dir: $outdir_full"
echo

## Download reference FASTA files, if necessary:
if [ ! -f "$refseq_dir"/pacific.01.pacific_9mers_nonGPU.h5 ]; then
    echo "Downloading PACIFIC reference data"
    wget -O "$refseq_dir"/pacific.01.pacific_9mers_nonGPU.h5 https://cloudstor.aarnet.edu.au/plus/s/Hwg20YRlua9a2OH/download
fi

## Run Pacific:
echo "## Running Pacific..."

python "$PACIFIC" \
  -i "$infile" \
  -m "$refseq_dir"/pacific.01.pacific_9mers_nonGPU.h5 \
  -t "$PACIFIC_DIR"/model/tokenizer.01.pacific_9mers.pickle \
  -l "$PACIFIC_DIR"/model/label_maker.01.pacific_9mers.pickle \
  -f fastq \
  -o "$outdir_full"

## Report:
echo "## Done with script."
date

