#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-bbmap-%j.out

## Report
echo "## Starting script bbmap.sh..."
date

## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/bbmap-env  # Also has samtools installed

## Bash strict mode
set -euo pipefail

## Command-line args
R1_in="$1"
R2_in="$2"
ref="$3"
bam_out="$4"

## Additional variables
n_cores="$SLURM_CPUS_ON_NODE"

## Process variables
bam_unsorted=${bam_out/.bam/}.unsorted.bam
outdir=$(dirname "$bam_out")
mkdir -p "$outdir"

## Report
echo "## Command-line args:"
echo "## Input FASTQ file - R1:          $R1_in"
echo "## Input FASTQ file - R2:          $R2_in"
echo "## Reference genome FASTQ:         $ref"
echo "## Output BAM file:                $bam_out"
echo -e "---------------------\n"

## Run BBMAP
echo "## Now running BBMap..."
bbmap.sh \
    ref="$ref" \
    in1="$R1_in" in2="$R2_in" \
    out="$bam_unsorted" \
    k=14 minid=0.9 build=1 \
    threads="$n_cores"

## Run samtools
echo -e "\n## Now running samtools sort..."
samtools sort "$bam_unsorted" > "$bam_out" && \
    rm "$bam_unsorted"

## Report
echo -e "\n## Output BAM file:"
ls -lh "$bam_out"
echo -e "\n## Done with script bbmap.sh"
date
echo
