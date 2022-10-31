#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=bowtie
#SBATCH --output=slurm-bowtie2-%j.out

## Load software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate bowtie2-env

## Bash strict mode
set -euo pipefail

## Command-line args
fa_in="$1"
fq_dir="$2"
outdir="$3"

## Check input
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments - you provided $#" && exit 1
[[ ! -d "$fq_dir" ]] && echo "## ERROR: Input FASTQ dir $fq_dir does not exist" && exit 1
[[ ! -f "$fa_in" ]] && echo "## ERROR: Input FASTA file $fa_in does not exist" && exit 1

## Process args
### Comma-delimited list of FASTQ files:
R1_list=$(echo "$fq_dir"/*R1*fastq.gz | sed 's/ /,/g')
R2_list=$(echo "$fq_dir"/*R2*fastq.gz | sed 's/ /,/g')

## Other parameters
n_cores=$SLURM_CPUS_PER_TASK
assembly_prefix=$(basename "$fa_in" |sed 's/.fa.*//')

## Report
echo
echo "## Starting script bowtie2.sh"
date
echo "## Args:                             $*"
echo
echo "## Input FASTQ dir:                  $fq_dir"
echo "## Input FASTA transcriptome file:   $fa_in"
echo "## Output dir:                       $outdir"
echo
echo "## List of R1 files:                 $R1_list"
echo "## List of R2 files:                 $R2_list"
echo -e "-------------------------------\n"

## Create output dir if needed
mkdir -p "$outdir"


# INDEX TRANSCRIPTOME ----------------------------------------------------------
echo -e "\n## Indexing the transcriptome..."
bowtie2-build \
    --threads "$n_cores" \
    "$fa_in" \
    $assembly_prefix

# MAP READS --------------------------------------------------------------------
echo -e "\n## Starting Bowtie2 mapping..."

bowtie2 -p 10 \
    -q \
    --no-unal \
    -k 20 \
    -x "$fa_in" \
    -1 "$R1_list" \
    -2 "$R2_list" \
    2>align_stats.txt |
    samtools view -@10 -Sb -o "$bam_out"


# WRAP-UP ----------------------------------------------------------------------
echo -e "\n-------------------------------"
echo "## Listing file in the output dir:"
ls -lh "$outdir"

echo -e "\n## Done with script trinity.sh"
date
