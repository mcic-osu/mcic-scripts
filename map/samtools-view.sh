#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --job-name=samtools-view
#SBATCH --output=slurm-samtools-view-%j.out

## Report
echo -e "\n## Starting script samtools-view.sh..."
date
echo

## Process args
bam_in=$1
outdir=$2
samtools_args=$3

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/samtools-env

## Define output file
bam_out="$outdir"/$(basename "$bam_in")

## Report
echo "## Input BAM file:                $bam_in"
echo "## Output BAM file:               $bam_out"
echo "## Other samtools arguments:      $samtools_args"
echo

## Make output dir
mkdir -p "$outdir"


# MAIN -------------------------------------------------------------------------
echo "## Now running samtools..."
samtools view "$bam_in" -b $samtools_args > "$bam_out"


# WRAP UP ----------------------------------------------------------------------
## Report
echo -e "\n## Output file:"
ls -lh "$bam_out"
echo -e "\n## Done with script samtools-view.sh"
date
echo
