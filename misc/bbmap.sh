#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-bbmap-%j.out

## Software
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
conda activate bbmap-env
conda activate --stack samtools-env

## Bash strict mode
set -euo pipefail

## Command-line args
R1_in="$1"
R2_in="$2"
scaffolds="$3"
bam_unsorted="$4"
bam_sorted="$5"

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"                            # Retrieve number of cores

## Other variables
outdir_unsorted=$(dirname "$bam_unsorted")
mkdir -p "$outdir_unsorted"
outdir_sorted=$(dirname "$bam_sorted")
mkdir -p "$outdir_sorted"

## Report
echo "## Starting script bbmap.sh..."
date
echo "## Command-line args:"
echo "## Input FASTQ file - R1: $R1_in"
echo "## Input FASTQ file - R2: $R2_in"
echo "## Assembled scaffold/contigs: $scaffolds"
echo "## Output BAM file - unsorted: $bam_unsorted"
echo "## Output BAM file - sorted: $bam_sorted"
echo -e "---------------------------\n\n"

## Run BBMAP
echo "Running BBMap..."
bbmap.sh \
    ref="$scaffolds" \
    in1="$R1_in" in2="$R2_in" \
    out="$bam_unsorted" \
    k=14 minid=0.9 build=1 \
    threads="$n_cores"

## Run samtools
echo -e "\nRunning samtools sort..."
samtools sort "$bam_unsorted" > "$bam_sorted"

## Report
echo -e "\n## Done with script bbmap.sh"
date
