#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --job-name=merge-bam
#SBATCH --output=slurm-merge-bam-%j.out

## Report
echo -e "\n## Starting script merge-bam.sh..."
date
echo

## Process args
indir=$1
bam_out=$2

## Load software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/samtools-env

## Determine file name for unsorted BAM file
bam_unsorted=${bam_out/.bam/_unsorted.bam}

## Report
echo "## Input dir:                     $indir"
echo "## Unsorted BAM file:             $bam_unsorted"
echo "## Output BAM file:               $bam_out"
echo

## Make output dir
outdir=$(dirname "$bam_out")
mkdir -p "$outdir"


# MAIN -------------------------------------------------------------------------
echo "## Now running samtools..."
samtools merge "$bam_unsorted" "$indir"/*bam
samtools sort "$bam_unsorted" > "$bam_out"

# WRAP UP ----------------------------------------------------------------------
## Report
echo -e "\n## Output file:"
ls -lh "$bam_out"
echo -e "\n## Done with script merge-bam.sh"
date
echo
