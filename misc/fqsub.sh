#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=12G
#SBATCH --time=30
#SBATCH --account=PAS0471
#SBATCH --job-name=fqsub
#SBATCH --output=slurm-fqsub-%j.out


# SET-UP & PARSE ARGS ----------------------------------------------------------
## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/seqtk-env

## Bash strict settings
set -euo pipefail

## Help function
Help()
{
   echo "## subsample_fastq.sh: script to subsample fastq files using seqtk"
   echo
   echo "## Syntax: subsample_fastq.sh -i <R1_in> -I <R2_in> -o <R1_out> -o <R2_out> [ -n <n_reads> | -p <prop_reads> ] [-h]"
   echo "## Options:"
   echo "## -h     Print help."
   echo "## -i     R1 input file (REQUIRED)"
   echo "## -I     R2 input file (REQUIRED)"
   echo "## -o     R1 output file (REQUIRED)"
   echo "## -O     R2 output file (REQUIRED)"
   echo "## -n     Number of reads (default: 100,000)"
   echo "## -p     Proportion of reads"
   echo
}

## Other parameters
random_seed=$RANDOM

## Option defaults
n_reads=100000
prop_reads="NA"

## Parse command-line options
while getopts ':i:I:o:O:n:p:h' flag; do
  case "${flag}" in
    i)  R1_in="$OPTARG" ;;
    I)  R2_in="$OPTARG" ;;
    o)  R1_out="$OPTARG" ;;
    O)  R2_out="$OPTARG" ;;
	n)  n_reads="$OPTARG" ;;
	p)  prop_reads="$OPTARG" ;;
    h)  Help && exit 0 ;;
	\?) echo "## ERROR: Invalid option" >&2 && exit 1 ;;
	:)  echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
  esac
done

## Check input - error out if neither n_reads or prop_reads is provided
[[ $n_reads = "NA" ]] && [[ $prop_reads = "NA" ]] && \
  echo "ERROR: neither a number of reads (-n) or a proportion of reads (-p) is provided" >&2 && exit 1
[[ ! -f "$R1_in" ]] && echo "ERROR: Input file $R1_in does not exist" >&2 && exit 1
[[ ! -f "$R2_in" ]] && echo "ERROR: Input file $R2_in does not exist" >&2 && exit 1

## Process parameters
### Number of reads in input FASTQ file
n_reads_total=$(zcat "$R1_in" | awk '{ s++ } END{ print s/4 }')
### If prop_reads is given, calculate n_reads
[[ $prop_reads != "NA" ]] && n_reads=$(python -c "print(int($n_reads_total * $prop_reads))")
### Create output dir if needed
outdir=$(dirname "$R1_in")
mkdir -p "$outdir"

## Report
echo -e "\n## Starting script subsample_fq.sh..."
date
echo "## Input R1:       $R1_in"
echo "## Input R2:       $R2_in"
echo "## Output R1:      $R1_out"
echo "## Output R2:      $R2_out"
echo "## Random seed:    $random_seed"
echo
echo "## Total (input) number of reads: $n_reads_total"
[[ $prop_reads != "NA" ]] && echo "## Proportion of reads to keep: $prop_reads"
echo "## Number of reads to keep: $n_reads"
echo
echo "## Listing input files:"
ls -lh "$R1_in"
ls -lh "$R2_in"
echo -e "----------------------------------\n\n"


# RUN SEQTK TO SUBSAMPLE FASTQ -------------------------------------------------
seqtk sample -s$random_seed "$R1_in" "$n_reads" | gzip > "$R1_out"
seqtk sample -s$random_seed "$R2_in" "$n_reads" | gzip > "$R2_out"


# REPORT AND FINALIZE --------------------------------------------------------
## Count reads in output
n_reads_R1_out=$(zcat "$R1_out" | awk '{ s++ } END{ print s/4 }')
n_reads_R2_out=$(zcat "$R2_out" | awk '{ s++ } END{ print s/4 }')
echo -e "\n## Number of reads in output files:"
echo "R1 out: $n_reads_R1_out"
echo "R2 out: $n_reads_R2_out"

## List output files
echo -e "\n----------------------------------"
echo -e "## Listing output files:"
ls -lh "$R1_out"
ls -lh "$R2_out"

echo -e "\n## Done with script."
date
