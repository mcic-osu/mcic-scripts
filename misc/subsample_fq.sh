#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10
#SBATCH --account=PAS0471

set -e -u -o pipefail

# Report:
echo -e "\n## Starting subsample_fastq script."
date

# SET-UP & PARSE ARGS ----------------------------------------------------------

# Help:
Help()
{
   # Display Help
   echo "## subsample_fastq.sh: script to subsample fastq files using seqtk"
   echo
   echo "## Syntax: subsample_fastq.sh -i <R1_in> -I <R2_in> -o <R1_out> -o <R2_out> [ -n <n_reads> | -p <prop_reads> ] [-h]"
   echo "## Options:"
   echo "## -h     Print help."
   echo "## -i     R1 input file (REQUIRED)"
   echo "## -I     R2 input file (REQUIRED)"
   echo "## -o     R1 output file (REQUIRED)"
   echo "## -O     R2 output file (REQUIRED)"
   echo "## -n     Number of reads"
   echo "## -p     Proportion of reads"
   echo
}

# Option defaults:
n_reads="NA"
prop_reads="NA"

# Get command-line options:
while getopts ':i:I:o:O:n:p:h' flag
do
  case "${flag}" in
    i)  R1_in="$OPTARG" ;;
    I)  R2_in="$OPTARG" ;;
    o)  R1_out="$OPTARG" ;;
    O)  R2_out="$OPTARG" ;;
	  n)  n_reads="$OPTARG" ;;
	  p)  prop_reads="$OPTARG" ;;
    h)  Help && exit 0 ;;
	  \?) echo "## trim.sh: ERROR: Invalid option" && exit 1 ;;
	  :)  echo "## trim.sh: ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

# Error out if neither n_reads or prop_reads is provided:
[[ $n_reads = "NA" ]] && [[ $prop_reads = "NA" ]] && \
  echo "Error: neither a number of reads (-n) or a proportion of reads (-p) is provided" && exit 1

# Number of reads in input fastq file:
n_reads_total=$(zcat $R1_in | awk '{ s++ } END{ print s/4 }')

# If prop_reads is given, calculate n_reads:
[[ $prop_reads != "NA" ]] && n_reads=$(python -c "print(int($n_reads_total * $prop_reads))")

# Random seed:
random_seed=$RANDOM

# REPORT -----------------------------------------------------------------------

echo "## Input R1: $R1_in"
echo "## Input R2: $R2_in"
echo "## Output R1: $R1_out"
echo "## Output R2: $R2_out"
echo
echo "## Total (input) number of reads: $n_reads_total"
[[ $prop_reads != "NA" ]] && echo "## Proportion of reads to keep: $prop_reads"
echo "## Number of reads to keep: $n_reads"
echo
echo "## Listing R1 input file:"
ls -lh $R1_in
echo "## Listing R2 input file:"
ls -lh $R2_in
echo
echo "## Random seed: $random_seed"
echo

# RUN SEQTK TO SUBSAMPLE FASTQ -------------------------------------------------

# Forward:
seqtk sample -s$random_seed $R1_in $n_reads | gzip > $R1_out

# Reverse:
seqtk sample -s$random_seed $R2_in $n_reads | gzip > $R2_out


# REPORT AND FINALIZE --------------------------------------------------------

# Count reads in output:
n_reads_R1_out=$(zcat $R1_out | awk '{ s++ } END{ print s/4 }')
n_reads_R2_out=$(zcat $R2_out | awk '{ s++ } END{ print s/4 }')
echo -e "\n## Number of reads in output files:"
echo "R1 out: $n_reads_R1_out"
echo "R2 out: $n_reads_R2_out"

# List files in output:
echo -e "\n## Listing output files:"
ls -lh $R1_out
ls -lh $R2_out

echo -e "\n## Done with script."
date
