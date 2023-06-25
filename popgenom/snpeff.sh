#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-snpeff-%j.out

## Load the software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/snpeff-env 

## Bash strict settings
set -ueo pipefail

## Process command line arguments
vcf_in=$1
vcf_out=$2
genomeID=$3

## Outdir
outdir=$(dirname "$vcf_out")
mkdir -p "$outdir"

## Report
echo "## Starting script snpeff.sh"
date
echo
echo "## Input VCF file:         $vcf_in"
echo "## Output VCF file:        $vcf_out"
echo "## Genome ID:              $genomeID"
echo -e "---------------------------\n"

## Make paths absolute and move into outdir
[[ ! $vcf_in =~ ^/ ]] && vcf_in="$PWD"/"$vcf_in"
[[ ! $vcf_out =~ ^/ ]] && vcf_out="$PWD"/"$vcf_out"
cd "$outdir" || exit

## Run SnpEff
snpEff ann \
    -o vcf \
    "$genomeID" \
    "$vcf_in" \
    >"$vcf_out"

## Report
echo "## Listing output files:"
ls -lh
echo
echo "## Done with script snpeff.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
