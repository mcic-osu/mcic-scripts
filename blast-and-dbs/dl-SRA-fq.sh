#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=120
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-download-fq-%j.out


# SETUP ------------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "## $0: Download FASTQ files from SRA."
    echo
    echo "## Syntax: $0 -i <SRA-ID> -I <SRA-ID-file> -o <output-dir>"
    echo 
    echo "## Required options:"
    echo "## -i STRING     A single SRA ID       [USE EITHER -i OR -I]"
    echo "## -I STRING     A file with SRA IDs   [USE EITHER -i OR -I]"
    echo "## -o STRING     Output dir"
    echo
    echo "## Other options:"
    echo "## -h     Print help this help message and exit"
    echo
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
    \?) echo "## $0: ERROR: Invalid option -$OPTARG" >&2 && exit 1 ;;
    :) echo "## $0: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Report
echo -e "\n## Starting script download-fq.sh..."
date
echo

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/sra-tools-env

## Bash strict settings
set -euo pipefail

## Check input
[[ "$SRA_ID" = "" ]] && [[ "$SRA_ID_file" = "" ]] && echo "ERROR: Must specify an SRA ID with -i or a file with SRA IDs with -I" >&2 && exit 1
[[ "$outdir" = "" ]] && echo "ERROR: must specify an output dir with -o" >&2 && exit 1

## Report
[[ "$SRA_ID" != "" ]] && echo "## SRA ID:               $SRA_ID"
echo "## File with SRA IDs:    $SRA_ID_file"
[[ "$SRA_ID_file" != "" ]] && echo "## Output dir:           $outdir"
echo

## Create output dir
mkdir -p "$outdir"


# DOWNLOAD FASTQ FILES ---------------------------------------------------------
echo "## Prefetching..."
prefetch "$SRA_ID" -O "$outdir"

echo "## Downloading..."
fasterq-dump "$SRA_ID" -O "$outdir"

echo "## Gzipping FASTQ files..."
find refdata/SRA -name "*${SRA_ID}*.fastq" -exec gzip {} \;


# WRAP UP ----------------------------------------------------------------------
echo -e "\n## Listing FASTQ files:"
ls -lh "$outdir"/*"${SRA_ID}"*fastq*
echo -e "\n## Done with script dl-SRA-fq.sh"
date
echo