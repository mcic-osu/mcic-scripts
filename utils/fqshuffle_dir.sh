#!/bin/bash


# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Constants
SUBSAMPLE_SCRIPT="mcic-scripts/misc/fqshuffle.sh"

## Help function
Help()
{
   echo "## fqshuffle_dir.sh: Script to shuffle FASTQ files using BBmap shuffle.sh"
   echo
   echo "## Syntax: fqshuffle_dir.sh -i <input-dir> -o <output-dir> [-h]"
   echo "## Options:"
   echo "## -h     Print help."
   echo "## -i     Input dir (REQUIRED)"
   echo "## -o     Output dir (REQUIRED)"
   echo "## -s     Sample ID pattern (to select only matching filenames)"
   echo
}

## Option defaults
sample_pattern=""

## Parse command-line options
while getopts ':i:o:s:h' flag; do
    case "${flag}" in
        i)  indir="$OPTARG" ;;
        o)  outdir="$OPTARG" ;;
        s)  sample_pattern="$OPTARG" ;;
        h)  Help && exit 0 ;;
	    \?) echo "## ERROR: Invalid option" >&2 && exit 1 ;;
	    :)  echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Check input
[[ ! -d "$indir" ]] && echo "## ERROR: Input dir $indir does not exist" >&2 && exit 1

## Report
echo
echo -e "\n## Starting script fqshuffle_dir.sh"
date
echo "## Input dir: $indir"
echo "## Output dir: $outdir"
echo "## Sample ID pattern: $sample_pattern"
echo -e "----------------------------------\n\n"

## Make output dir if needed
mkdir -p "$outdir"


# RUN SCRIPT FOR EACH PAIR OF FASTQ FILES --------------------------------------
for R1 in "$indir"/${sample_pattern}*_R1*.fastq.gz; do
    R1=$(basename "$R1")
    R2=${R1/_R1/_R2}

    ## Report input files:
    echo "## R1 input file:"
    ls -lh "$indir/$R1"
    echo "## R2 input file:"
    ls -lh "$indir/$R2"

    ## SLURM log file
    sample_id=$(basename "$R1" .fastq.gz)
    log=slurm-fqshuffle-"$sample_id"-%j.out

    sbatch -o "$log" "$SUBSAMPLE_SCRIPT" \
        -i "$indir/$R1" -I "$indir/$R2" \
        -o "$outdir/$R1" -O "$outdir/$R2"
  
    echo -e "---------------------------\n"
done


# REPORT AND FINALIZE --------------------------------------------------------
echo -e "## Done with script fqshuffle_dir.sh"
date