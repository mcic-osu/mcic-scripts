#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-trim-%j.out

# Software:
source ~/.bashrc
source activate trimmomatic-env

set -euo pipefail

n_cores="$SLURM_CPUS_ON_NODE"

# Help:
Help() {
    # Display Help
    echo "## script to trim sequences in FASTQC files using Trimmomatic"
    echo
    echo "## Syntax: trim.sh -b barcode-file -i input-dir -l library-ID -o output-dir -s steps-to-skip [-h]"
    echo "## Options:"
    echo "## -h     Print help."
    echo "## -i     Input R1 sequence file (REQUIRED)"
    echo "## -o     Output dir (default: results/trim/)"
    echo "## -t     Trimming-stats dir (default: $outdir/log)"
    echo "## -a     Adapter file (default: none. Other options: NexteraPE-PE.fa, TruSeq2-PE.fa, TruSeq3-PE.fa)"
    echo
}

## Option defaults:
R1_in=""
outdir="results/trim"
stats_dir="$outdir/log"
adapter_file="NA"

## Parse options:
while getopts ':i:o:t:a:h' flag; do
    case "${flag}" in
    i) R1_in="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    t) stats_dir="$OPTARG" ;;
    a) adapter_file="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option" && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

[[ -z "$R1_in" ]] && echo "## ERROR: Please provide R1 input FASTQ file with -i flag" && exit 1

## Process args:
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?[0-9]).*fastq.gz/\1/')

R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}

R1_basename=$(basename "$R1_in" .fastq.gz)
R2_basename=$(basename "$R2_in" .fastq.gz)

sample_ID=${R1_basename/"$R1_suffix"/}

discard_dir="$outdir"/discard                   # Dir for discarded sequences
trimstats_file="$stats_dir"/"$sample_ID".trimstats.txt # File with Trimmomatic stdout

R1_out="$outdir"/"$R1_basename".fastq.gz # Output R1 file
R2_out="$outdir"/"$R2_basename".fastq.gz # Output R2 file

R1_discard=$discard_dir/"$R1_basename"_U1.fastq.gz # Output file for discarded R1 reads
R2_discard=$discard_dir/"$R2_basename"_U2.fastq.gz # Output file for discarded R2 reads

## Adapter argument:
if [[ $adapter_file = "NA" ]]; then
    adapter_arg=""
else
    adapter_arg=" ILLUMINACLIP:$adapter_file:2:30:10"
fi

## Report:
echo "## Starting script trimmomatic.sh"
date
echo "## Input R1: $R1_in"
echo "## Input R2: $R2_in"
echo "## Output R1: $R1_out"
echo "## Output R2: $R2_out"
echo "## File for discarded R1 seqs: $R1_discard"
echo "## File for discarded R2 seqs: $R2_discard"
echo "## Adapter argument: $adapter_arg"
echo

[[ ! -d $discard_dir ]] && mkdir -p "$discard_dir" # Create dir for discarded sequences if it doesn't exist
[[ ! -d $stats_dir ]] && mkdir -p "$stats_dir"     # Create dir for stdout file if it doesn't exist

[[ ! -f $R1_in ]] && echo "## ERROR: Input file R1_in does not exist" && exit 1
[[ ! -f $R2_in ]] && echo "## ERROR: Input file R2_in does not exist" && exit 1

echo "## Listing input files:"
ls -lh "$R1_in"
ls -lh "$R2_in"
echo

## Run Trimmomatic:
trimmomatic PE -threads "$n_cores" -phred33 \
    "$R1_in" "$R2_in" "$R1_out" "$R1_discard" "$R2_out" "$R2_discard"$adapter_arg \
    AVGQUAL:28 LEADING:20 TRAILING:20 MINLEN:36 \
    2>&1 | tee "$trimstats_file"

## Report:
nreads_raw=$(zcat "$R1_in" | awk '{ s++ } END{ print s/4 }')
nreads_trim=$(zcat "$R1_out" | awk '{ s++ } END{ print s/4 }')
echo "Number of raw / trimmed read-pairs: $nreads_raw / $nreads_trim"

echo -e "\n## Listing output files:"
ls -lh "$R1_out"
ls -lh "$R2_out"

echo "## Done with script."
date
