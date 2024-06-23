#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-bam_idx-%j.out

## Report
echo -e "\n## Starting script bam_idx.sh..."
date
echo

## Process args
bam=$1

## Load software
module load python/3.6-conda5.2
source activate /fs/ess/PAS0471/jelmer/conda/samtools

## Report
echo "## Input BAM file:      $bam"
echo


# MAIN -------------------------------------------------------------------------
samtools index "$bam"


# WRAP UP ----------------------------------------------------------------------
## Report
echo -e "\n## Output files:"
ls -lh "$bam"*
echo -e "\n## Done with script idx_bam.sh"
date
echo
