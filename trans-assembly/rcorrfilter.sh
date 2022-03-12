#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=readfilter
#SBATCH --output=slurm-readfilter-%j.out

## Software and scripts
source ~/.bashrc
[[ $(which conda) = ~/miniconda3/bin/conda ]] || module load python/3.6-conda5.2
source activate rcorrector-env

##> Script from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
script=software/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py

## Bash strict mode
set -euo pipefail

## Command-line args
R1="$1"
outdir="$2"

## Process args
[[ "$#" != 2 ]] && echo "## ERROR: Please provide 2 arguments - you provided $#" && exit 1
R2=${R1/_R1/_R2}
sample_id=$(echo "$(basename $R1)" | sed 's/_R1.*//')

R1_out="$outdir"/$(basename "$R1" .cor.fq.gz).fastq
R2_out="$outdir"/$(basename "$R2" .cor.fq.gz).fastq

## Report
echo
echo "## Starting script discard_bad_reads.sh"
date
echo "## R1 input file:             $R1"
echo "## R2 input file (inferred):  $R2"
echo "## Output dir:                $outdir"
echo "## R1 output file:            $R1_out"
echo "## R2 output file:            $R2_out"
echo "## Sample ID (inferred):      $sample_id"
echo -e "-------------------------------\n"

## Create output dirs if needed
mkdir -p "$outdir"


# RUN PYTHON SCRIPT ------------------------------------------------------------
echo "## Running read-filter script..."
python "$script" -s "$sample_id" -1 "$R1" -2 "$R2"


# WRAP-UP ----------------------------------------------------------------------
## Move output files
echo -e "\n## Moving and gzipping output files..."
mv unfixrm*"$sample_id"*R1*cor.fq "$R1_out"
gzip -f "$R1_out"
mv unfixrm*"$sample_id"*R2*cor.fq "$R2_out"
gzip -f "$R2_out"

mv rmunfixable_"$sample_id".log "$outdir"

echo -e "\n-------------------------------"
echo "## Listing output files:"
ls -lh "$R1_out".gz "$R2_out".gz

echo -e "\n## Done with script readfilter.sh"
date



