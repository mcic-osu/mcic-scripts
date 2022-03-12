#!/bin/bash

# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Constants
SUBSAMPLE_SCRIPT="mcic-scripts/misc/fqsub.sh"

## Help function
Help()
{
   echo "## fqsub_dir.sh: script to subsample a dir of FASTQ files using seqtk"
   echo
   echo "## Syntax: fqsub_dir.sh -i <input-dir> -o <output-dir> ..."
   echo
   echo "## Required options:"
   echo "## -i     Input dir with FASTQ files"
   echo "## -o     Output dir for subsamples FASTQ files"
   echo
   echo "## Other options:"
   echo "## -n     Number of reads (default: 100,000)"
   echo "## -p     Proportion of reads (alternative to -n)"
   echo "## -s     Sample ID pattern (to select only matching filenames)"
   echo "## -h     Print help."
   echo
}

## Option defaults
n_reads=100000
prop_reads="NA"
sample_pattern=""

## Parse command-line options
while getopts ':i:o:n:p:s:h' flag; do
    case "${flag}" in
        i)  indir="$OPTARG" ;;
        o)  outdir="$OPTARG" ;;
	    n)  n_reads="$OPTARG" ;;
	    p)  prop_reads="$OPTARG" ;;
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
echo -e "\n## Starting script fqsub_dir.sh"
date
echo "## Input dir:                    $indir"
echo "## Output dir:                   $outdir"
echo "## Number of reads to keep:      $n_reads"
echo "## Proportion of reads to keep:  $prop_reads"
echo "## Sample ID pattern:            $sample_pattern"
echo -e "------------------------\n"

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
    log=slurm-fqsub-"$sample_id"-%j.out

    sbatch -o "$log" "$SUBSAMPLE_SCRIPT" \
        -i "$indir/$R1" -I "$indir/$R2" \
        -o "$outdir/$R1" -O "$outdir/$R2" \
        -n "$n_reads" -p "$prop_reads"
  
    echo -e "---------------------------\n"
done


# REPORT AND FINALIZE --------------------------------------------------------
echo -e "## Done with script fqsub_dir.sh"
date