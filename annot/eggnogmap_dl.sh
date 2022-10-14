#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --job-name=eggnogmap_dl
#SBATCH --output=slurm-eggnogmap-dl-%j.out

## Arguments
outdir=$1

## Software
module load python/3.6-conda5.2
source activate /users/PAS0471/jelmer/miniconda3/envs/eggnogg-env

## Bash strict mode
set -euo pipefail

## Make output dir
mkdir -p "$outdir"

## Report
echo "## Starting script eggnogmap_dl.sh"
date
echo
echo "## Output dir:                           $outdir"
echo -e "--------------------\n"

## Download
download_eggnog_data.py \
    -y \
    -P \
    -M \
    --data_dir "$outdir"


## Report
echo -e "\n-------------------------------"
echo "## Listing files in the output dir:"
ls -lh "$outdir"
echo -e "\n## Done with script eggnogmap_dl.sh"
date
echo
sacct -j "$SLURM_JOB_ID" -o JobID,AllocTRES%50,Elapsed,CPUTime,TresUsageInTot,MaxRSS
echo
