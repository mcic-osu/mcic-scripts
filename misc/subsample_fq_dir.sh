#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120
#SBATCH --account=PAS0471


# SETUP ------------------------------------------------------------------------
## Bash strict settings
set -euo pipefail

## Constants
SUBSAMPLE_SCRIPT="mcic-scripts/misc/subsample_fq.sh"

## Help function
Help()
{
   # Display Help
   echo "## subsample_fastq.sh: script to subsample fastq files using seqtk"
   echo
   echo "## Syntax: subsample_fastq.sh -i <R1_in> -I <R2_in> -o <R1_out> -o <R2_out> [ -n <n_reads> | -p <prop_reads> ] [-h]"
   echo "## Options:"
   echo "## -h     Print help."
   echo "## -i     Input dir (REQUIRED)"
   echo "## -o     Output dir (REQUIRED)"
   echo "## -n     Number of reads"
   echo "## -p     Proportion of reads"
   echo
}

## Option defaults
n_reads="NA"
prop_reads="NA"

## Parse command-line options
while getopts ':i:o:n:p:h' flag
do
  case "${flag}" in
    i)  indir="$OPTARG" ;;
    o)  outdir="$OPTARG" ;;
	  n)  n_reads="$OPTARG" ;;
	  p)  prop_reads="$OPTARG" ;;
    h)  Help && exit 0 ;;
	  \?) echo "## trim.sh: ERROR: Invalid option" && exit 1 ;;
	  :)  echo "## trim.sh: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Report
echo -e "\n## Starting script subsample_fastq_dir.sh"
date
echo "## Input dir: $indir"
echo "## Output dir: $outdir"
echo "## Number of reads to keep: $n_reads"
echo "## Proportion of reads to keep: $prop_reads"
echo -e "----------------------------------\n\n"

## Make output dir if needed
mkdir -p "$outdir"


# RUN SCRIPT FOR EACH PAIR OF FASTQ FILES --------------------------------------
for R1 in "$indir"/*_R1*.fastq.gz; do
    R1=$(basename "$R1")
    R2=${R1/_R1/_R2}

    ## Report input files:
    echo "## R1 input file:"
    ls -lh "$indir/$R1"
    echo "## R2 input file:"
    ls -lh "$indir/$R2"

    ## SLURM log file
    sample_id=$(basename "$R1" .fastq.gz)
    log=slurm-subsample-fastq_"$sample_id"_%j.out

    sbatch -o "$log" "$SUBSAMPLE_SCRIPT" \
        -i "$indir/$R1" -I "$indir/$R2" \
        -o "$outdir/$R1" -O "$outdir/$R2" \
        -n "$n_reads" -p "$prop_reads"
  
    echo -e "\n--------------------------------------------------------\n"
done


# REPORT AND FINALIZE --------------------------------------------------------
echo -e "\n\n## Done with script subsample_fastq_dir.sh."
date