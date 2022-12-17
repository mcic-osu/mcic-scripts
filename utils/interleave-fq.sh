#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --output=slurm-interleave-%j.out

## Parse positional args
R1_in=$1
outdir=$2

## Software
module load python
source activate /fs/ess/PAS0471/jelmer/conda/bbmap-38.96

## Bash strict settings
set -euo pipefail

## Process parameters
R1_suffix=$(echo "$R1_in" | sed -E 's/.*(_R?1).*fa?s?t?q.gz/\1/')
R2_suffix=${R1_suffix/1/2}
R2_in=${R1_in/$R1_suffix/$R2_suffix}
R1_basename=$(basename "$R1_in" .fastq.gz)
sampleID=${R1_basename/"$R1_suffix"/}
out="$outdir"/"$sampleID".fastq.gz

## Report
echo "## Starting script interleave.sh"
date
echo
echo "## R1 input file:              $R1_in"
echo "## R2 input file:              $R2_in"
echo "## Interleave output file:     $out"
echo

echo "## Now running reformat.sh..."
reformat.sh in1="$R1_in" in2="$R2_in" out="$out"

echo -e "\n## Done with script interleave.sh"
date
echo
