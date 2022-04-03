#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm-download-refseq-%j.out

## Bash strict settings
set -euo pipefail

## Parameters
refseq_lib=$1     # Options are the dir names here: https://ftp.ncbi.nlm.nih.gov/refseq/release/
outdir=$2

## Report
echo -e "\n## Starting script dl-refseq.sh..."
date
echo
echo "## Refseq library:         $refseq_lib"
echo "## Output dir:             $outdir"
echo

## NCBI key
export NCBI_API_KEY=34618c91021ccd7f17429b650a087b585f08
#? https://ftp.ncbi.nlm.nih.gov/refseq/release/
#? See also https://github.com/DerrickWood/kraken2/issues/115

## Make outdir
mkdir -p "$outdir"

## Download
wget \
    ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"$refseq_lib"/*genomic.fna.gz \
    -P "$outdir"

## Report
echo -e "\n## Listing files in output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script dl-refseq.sh..."
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
