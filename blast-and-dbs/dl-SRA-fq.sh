#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=120
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-download-fq-%j.out

# Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate sra-tools-env

## Bash strict settings
set -euo pipefail

## Help
Help() {
    # Display Help
    echo
    echo "## $0: Download FASTQ files from SRA."
    echo
    echo "## Syntax: $0 -i <input-sequence-file> -I <> -o <output-dir> [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -i     A single SRA ID"
    echo "## -I     A file with SRA IDs"
    echo "## -o     Output dir (REQUIRED)"
    echo "## Example: $0 -i SRR3678332 -o refdata/SRA"
    echo "## To submit the OSC queue, preface with 'sbatch': sbatch $0 ..."
    echo
}

## Option defaults
SRA_ID=""
SRA_ID_file=""
outdir=""

## Get command-line options
while getopts 'i:I:o:h' flag; do
    case "${flag}" in
    i) SRA_ID="$OPTARG" ;;
    I) SRA_ID_file="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## $0: ERROR: Invalid option" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Process options
[[ "$SRA_ID" = "" ]] && [[ "$SRA_ID_file" = "" ]] && echo "ERROR: Must specify an SRA ID with -i or a file with SRA IDs with -I" && exit 1
[[ "$outdir" = "" ]] && echo "ERROR: must specify an output dir with -o" && exit 1

## Report
echo "## Starting script download-fq.sh..."
date
echo "## Command-line args:"
[[ "$SRA_ID" != "" ]] && echo "## SRA ID: $SRA_ID"
[[ "$SRA_ID_file" != "" ]] && echo "## File with SRA IDs: $SRA_ID_file"
echo "## Output dir: $outdir"
echo

## Create output dir
mkdir -p "$outdir"


# DOWNLOAD FASTQ FILES ---------------------------------------------------------
## Download FASTQ files
echo "## Prefetching..."
prefetch "$SRA_ID" -O "$outdir"

echo "## Downloading..."
fasterq-dump "$SRA_ID" -O "$outdir"

## gzip FASTQ files
echo "## Gzipping FASTQ files..."
find refdata/SRA -name "*${SRA_ID}*.fastq" -exec gzip {} \;

## Report
echo "## Listing FASTQ files:"
ls -lh "$outdir"/*"${SRA_ID}"*fastq*

echo -e "\n## Done with script."
date