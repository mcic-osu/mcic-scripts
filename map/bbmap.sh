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
scaffolds="$2"
bam_sorted="$3"

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"                            # Retrieve number of cores

R2_in=${R1_in/_R1/_R2}
bam_unsorted=${bam_sorted/.bam/.unsorted.bam}

outdir=$(dirname "$bam_sorted")
mkdir -p "$outdir"

## Report
echo "## Starting script bbmap.sh..."
date
echo "## Command-line args:"
echo "## Input FASTQ file - R1: $R1_in"
echo "## Assembled scaffold/contigs: $scaffolds"
echo "## Output BAM file - sorted: $bam_sorted"
echo
echo "## Other (inferred) variables:"
echo "## Input FASTQ file - R2: $R2_in"
echo "## Output BAM file - unsorted: $bam_unsorted"
echo -e "---------------------------\n\n"

## Test
[[ ! -f "$R1_in" ]] && echo "## ERROR: Input file R1_in $R1_in does not exist" && exit 1
[[ ! -f "$R2_in" ]] && echo "## ERROR: Input file R2_in $R2_in does not exist" && exit 1

## Run BBMAP
echo "Running BBMap..."
bbmap.sh \
    ref="$scaffolds" \
    in1="$R1_in" in2="$R2_in" \
    out="$bam_unsorted" \
    k=14 minid=0.9 build=1 \
    threads="$n_cores"

## Run samtools
echo -e "\n## Running samtools sort..."
samtools sort "$bam_unsorted" > "$bam_sorted"

## Report
echo -e "\n## Listing output BAM files:"
ls -lh "$bam_unsorted"
ls -lh "$bam_sorted"

echo -e "\n## Done with script bbmap.sh"
date
