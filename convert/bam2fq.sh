#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-bam2fq-%j.out

## Report
echo -e "\n## Starting script bam2fq.sh..."
date
echo

## Process args
bam=$1
outdir=$2

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bedtools-env

## Output files
mkdir -p "$outdir"
sampleID=$(basename "$bam" .bam)
fq_R1=$outdir/"$sampleID"_R1.fastq
fq_R2=$outdir/"$sampleID"_R2.fastq

## Report
echo "## Input BAM file:      $bam"
echo "## Output R1 FASTQ:     $fq_R1"
echo "## Output R2 FASTQ:     $fq_R2"
echo

# MAIN -------------------------------------------------------------------------
## Convert BAM to FASTQ
bedtools bamtofastq -i "$bam" -fq "$fq_R1" -fq2 "$fq_R2"

## Gzip FASTQ files
gzip "$fq_R1"
gzip "$fq_R2"


# WRAP UP ----------------------------------------------------------------------
## Report
echo -e "\n## Output files:"
ls -lh "$fq_R1".gz "$fq_R2".gz
echo -e "\n## Done with script bam2fq.sh"
date
echo
