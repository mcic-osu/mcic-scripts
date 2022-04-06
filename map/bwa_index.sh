#!/bin/bash

#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm-bwa-%j.out

## Report
echo -e "\n## Starting script bwa_index.sh"
date
echo

## Load software
module load python/3.6-conda5.2
source activate /fs/project/PAS0471/jelmer/conda/bwa-0.7.17

## Bash strict settings
set -euo pipefail

## Command-line args
ref=$1

## Index the FASTA file
bwa index "$ref"

## Report
echo -e "\n## Output dir:"
ls -lh "$(dirname "$ref")"
echo -e "\n## Done with script bwa_index.sh"
date
echo
