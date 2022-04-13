#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm-trimmomatic-%j.out


# SETUP ------------------------------------------------------------------------
## Help function
Help() {
    echo
    echo "## $0: Script to trim sequences in FASTQ files using Trimmomatic"
    echo
    echo "## Syntax: $0 -i <input R1 FASTQ> [ -o <output dir> ] [ -a <adapter-file ] [ -p <trim-params> ] [-h]"
    echo "## Options:"
    echo "## -h         Print help."
    echo "## -i STR     Input R1 sequence file (REQUIRED)"
    echo "## -o STR     Output dir (default: 'results/trimmomatic/')"
    echo "## -a STR     Adapter file (default: none. Other options: NexteraPE-PE.fa, TruSeq2-PE.fa, TruSeq3-PE.fa)"
    echo "              Or provide your own FASTA file with adapters, e.g. 'adapters.fa' from BBduk"
    echo "## -p STR     Trimming paramaters for Trimmomatic"
    echo "              (default: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36')"
    echo
    echo "## Example command:"
    echo "$0 -i data/fastq/A1_L001_R1_001.fastq.gz -o results/trimmomatic"
    echo
}

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/trimmomatic-env

## Bash strict settings
set -euo pipefail

## Option defaults
R1_in=""
outdir="results/trimmomatic"
adapter_file="NA"

trim_param="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
# Same as example on https://github.com/usadellab/Trimmomatic and https://rpubs.com/ednachiang/MetaG_Pipeline
# Alternatively, example of a much stricter mode: "AVGQUAL:28 LEADING:20 TRAILING:20 MINLEN:36"

## Parse options
while getopts ':i:o:a:p:h' flag; do
    case "${flag}" in
    i) R1_in="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    a) adapter_file="$OPTARG" ;;
    p) trim_param="$OPTARG" ;;
    h) Help && exit 0 ;;
    \?) echo "## ERROR: Invalid option" && exit 1 ;;
    :) echo "## ERROR: Option -$OPTARG requires an argument." >&2 && exit 1 ;;
    esac
done

## Other parameters
n_cores="$SLURM_CPUS_ON_NODE"

## Process parameters
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?[0-9]).*fastq.gz/\1/')
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}

R1_basename=$(basename "$R1_in" .fastq.gz)
R2_basename=$(basename "$R2_in" .fastq.gz)

sample_ID=${R1_basename/"$R1_suffix"/}

stats_dir="$outdir/log"
discard_dir="$outdir"/discard                          # Dir for discarded sequences
trimstats_file="$stats_dir"/"$sample_ID".trimstats.txt # File with Trimmomatic stdout

R1_out="$outdir"/"$R1_basename".fastq.gz # Output R1 file
R2_out="$outdir"/"$R2_basename".fastq.gz # Output R2 file

R1_discard=$discard_dir/"$R1_basename"_U1.fastq.gz # Output file for discarded R1 reads
R2_discard=$discard_dir/"$R2_basename"_U2.fastq.gz # Output file for discarded R2 reads

## Adapter argument
if [[ $adapter_file = "NA" ]]; then
    adapter_arg=""
else
    adapter_arg=" ILLUMINACLIP:$adapter_file:2:30:10:2:True"
    # As in the example here https://github.com/usadellab/Trimmomatic
fi

## Trimming parameters arg
trim_param_arg=" $trim_param"

## Report
echo "## Starting script trimmomatic.sh"
date
echo
echo "## Input R1:                     $R1_in"
echo "## Output dir:                   $outdir"
echo "## Trimming params:              $trim_param"
echo "## Adapter argument:             $adapter_arg"
echo
echo "## Input R2:                     $R2_in"
echo "## Output R1:                    $R1_out"
echo "## Output R2:                    $R2_out"
echo "## File for discarded R1 seqs:   $R1_discard"
echo "## File for discarded R2 seqs:   $R2_discard"
echo
echo "## Listing input files:"
ls -lh "$R1_in"
ls -lh "$R2_in"
echo -e "--------------------------------------\n"

## Check input
[[ -z "$R1_in" ]] && echo "## ERROR: Please provide R1 input FASTQ file with -i flag" && exit 1
[[ ! -f $R1_in ]] && echo "## ERROR: Input file R1_in ($R1_in) does not exist" && exit 1
[[ ! -f $R2_in ]] && echo "## ERROR: Input file R2_in ($R2_in) does not exist" && exit 1

## Create output dirs
mkdir -p "$discard_dir"   # Create dir for discarded sequences if it doesn't exist
mkdir -p "$stats_dir"     # Create dir for stdout file if it doesn't exist


# RUN TRIMMOMATIC --------------------------------------------------------------
echo -e "## Starting Trimmomatic run..."
trimmomatic PE \
    -threads "$n_cores" \
    "$R1_in" "$R2_in" \
    "$R1_out" "$R1_discard" \
    "$R2_out" "$R2_discard"${adapter_arg}${trim_param_arg} \
    2>&1 | tee "$trimstats_file"


# WRAP UP ----------------------------------------------------------------------
nreads_raw=$(zcat "$R1_in" | awk '{ s++ } END{ print s/4 }')
nreads_trim=$(zcat "$R1_out" | awk '{ s++ } END{ print s/4 }')
echo -e "\n## Number of raw / trimmed read-pairs: $nreads_raw / $nreads_trim"

echo -e "\n## Listing output files:"
ls -lh "$R1_out"
ls -lh "$R2_out"

echo "## Done with script trimmomatic.sh"
date
